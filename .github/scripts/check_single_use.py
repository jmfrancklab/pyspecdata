#!/usr/bin/env python3
# In GitHub Actions, this checker is intended to run only on changed Python
# files, not on the whole repository.
#
# Local usage before pushing a branch:
#   python .github/scripts/check_single_use.py path/to/file.py
#
# You can also pass several files explicitly:
#   python .github/scripts/check_single_use.py path/to/one.py path/to/two.py
#
# This is useful when you want to check a file of your choice before opening a
# PR.
"""Check changed files for single-use module variables and single-call
functions."""

from __future__ import annotations

import argparse
import ast
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re
import sys


@dataclass(frozen=True)
class Violation:
    """A single finding from the checker."""

    path: str
    line: int
    name: str
    kind: str
    definition_line: int
    use_line: int


class ModuleAnalyzer(ast.NodeVisitor):
    """Collect module-level definitions and direct calls."""

    def __init__(self):
        self.scope_depth = 0
        self.class_depth = 0
        self.main_guard_depth = 0
        self.module_variable_definitions = defaultdict(list)
        self.for_target_definitions = set()
        self.multi_name_lhs_definitions = set()
        self.with_alias_definitions = set()
        self.module_variable_uses = defaultdict(list)
        self.top_level_function_definitions = defaultdict(list)
        self.direct_function_calls = defaultdict(list)
        self.main_guard_direct_calls = defaultdict(list)

    def _record_target(
        self,
        node: ast.AST,
        allow_for_target_exception: bool = False,
        allow_multi_name_lhs_exception: bool = False,
        allow_with_alias_exception: bool = False,
    ):
        target_names = []
        nodes_to_visit = [node]
        while nodes_to_visit:
            current_node = nodes_to_visit.pop()
            if isinstance(current_node, ast.Name):
                target_names.append((current_node.id, current_node.lineno))
                continue
            if isinstance(current_node, (ast.Tuple, ast.List)):
                nodes_to_visit.extend(reversed(current_node.elts))
                continue
            if isinstance(current_node, ast.Starred):
                nodes_to_visit.append(current_node.value)
        for name, lineno in target_names:
            self.module_variable_definitions[name].append(lineno)
        if allow_for_target_exception:
            self.for_target_definitions.update(target_names)
        if (
            allow_multi_name_lhs_exception
            and isinstance(node, (ast.Tuple, ast.List))
            and len(target_names) > 1
        ):
            self.multi_name_lhs_definitions.update(target_names)
        if allow_with_alias_exception:
            self.with_alias_definitions.update(target_names)

    def _visit_nested_scope(self, node):
        self.scope_depth += 1
        self.generic_visit(node)
        self.scope_depth -= 1

    def visit_Name(self, node: ast.Name):
        if (
            self.scope_depth == 0
            and self.class_depth == 0
            and isinstance(node.ctx, ast.Load)
        ):
            self.module_variable_uses[node.id].append(node.lineno)

    def visit_Assign(self, node: ast.Assign):
        if self.scope_depth == 0 and self.class_depth == 0:
            for target in node.targets:
                self._record_target(
                    target, allow_multi_name_lhs_exception=True
                )
        self.generic_visit(node)

    def visit_AnnAssign(self, node: ast.AnnAssign):
        if self.scope_depth == 0 and self.class_depth == 0:
            self._record_target(node.target)
        self.generic_visit(node)

    def visit_AugAssign(self, node: ast.AugAssign):
        if self.scope_depth == 0 and self.class_depth == 0:
            self._record_target(node.target)
        self.generic_visit(node)

    def visit_For(self, node: ast.For):
        if self.scope_depth == 0 and self.class_depth == 0:
            self._record_target(node.target, allow_for_target_exception=True)
        self.generic_visit(node)

    def visit_AsyncFor(self, node: ast.AsyncFor):
        if self.scope_depth == 0 and self.class_depth == 0:
            self._record_target(node.target, allow_for_target_exception=True)
        self.generic_visit(node)

    def visit_With(self, node: ast.With):
        if self.scope_depth == 0 and self.class_depth == 0:
            for item in node.items:
                if item.optional_vars is not None:
                    self._record_target(
                        item.optional_vars, allow_with_alias_exception=True
                    )
        self.generic_visit(node)

    def visit_AsyncWith(self, node: ast.AsyncWith):
        if self.scope_depth == 0 and self.class_depth == 0:
            for item in node.items:
                if item.optional_vars is not None:
                    self._record_target(
                        item.optional_vars, allow_with_alias_exception=True
                    )
        self.generic_visit(node)

    def visit_NamedExpr(self, node: ast.NamedExpr):
        if self.scope_depth == 0 and self.class_depth == 0:
            self._record_target(node.target)
        self.generic_visit(node)

    def visit_ExceptHandler(self, node: ast.ExceptHandler):
        if (
            self.scope_depth == 0
            and self.class_depth == 0
            and node.name is not None
        ):
            self.module_variable_definitions[node.name].append(node.lineno)
        self.generic_visit(node)

    def visit_ClassDef(self, node: ast.ClassDef):
        self.class_depth += 1
        self.generic_visit(node)
        self.class_depth -= 1

    def visit_If(self, node: ast.If):
        in_main_guard = (
            self.scope_depth == 0
            and self.class_depth == 0
            and isinstance(node.test, ast.Compare)
            and isinstance(node.test.left, ast.Name)
            and node.test.left.id == "__name__"
            and len(node.test.ops) == 1
            and isinstance(node.test.ops[0], ast.Eq)
            and len(node.test.comparators) == 1
            and isinstance(node.test.comparators[0], ast.Constant)
            and node.test.comparators[0].value == "__main__"
        )
        if in_main_guard:
            self.main_guard_depth += 1
        self.generic_visit(node)
        if in_main_guard:
            self.main_guard_depth -= 1

    def visit_FunctionDef(self, node: ast.FunctionDef):
        if self.scope_depth == 0 and self.class_depth == 0:
            self.top_level_function_definitions[node.name].append(node.lineno)
        self._visit_nested_scope(node)

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef):
        if self.scope_depth == 0 and self.class_depth == 0:
            self.top_level_function_definitions[node.name].append(node.lineno)
        self._visit_nested_scope(node)

    def visit_Lambda(self, node: ast.Lambda):
        self._visit_nested_scope(node)

    def visit_Call(self, node: ast.Call):
        if isinstance(node.func, ast.Name):
            self.direct_function_calls[node.func.id].append(node.lineno)
            if self.main_guard_depth > 0 and self.scope_depth == 0:
                self.main_guard_direct_calls[node.func.id].append(node.lineno)
        self.generic_visit(node)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Check for module-level variables that are defined once and used "
            "once, and top-level functions that are directly called once."
        )
    )
    parser.add_argument("paths", nargs="+", help="Python files to inspect")

    fold_start = re.compile(r"^\s*#\s*\{\{\{")
    fold_end = re.compile(r"^\s*#\s*\}\}\}")
    violations = []
    parse_failures = []

    for raw_path in parser.parse_args(argv).paths:
        path = Path(raw_path)
        try:
            source = path.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=str(path))
        except (OSError, SyntaxError, UnicodeDecodeError) as exc:
            parse_failures.append(
                f"{path}: could not analyze this file: {exc}"
            )
            continue

        lines = source.splitlines()
        preamble_limit = len(lines)
        for stmt in tree.body:
            if (
                isinstance(stmt, ast.Expr)
                and isinstance(stmt.value, ast.Constant)
                and isinstance(stmt.value.value, str)
            ):
                continue
            if isinstance(
                stmt, (ast.Import, ast.ImportFrom, ast.Assign, ast.AnnAssign)
            ):
                continue
            preamble_limit = stmt.lineno - 1
            break

        allowed_blocks = []
        line_idx = 0
        while line_idx < len(lines):
            if not fold_start.match(lines[line_idx]):
                line_idx += 1
                continue
            block_start = line_idx + 1
            block_lines = [lines[line_idx]]
            scan_idx = line_idx + 1
            while scan_idx < len(lines):
                block_lines.append(lines[scan_idx])
                if fold_end.match(lines[scan_idx]):
                    block_end = scan_idx + 1
                    if block_start <= preamble_limit:
                        opening_comment_lines = []
                        for line in block_lines:
                            stripped = line.strip()
                            if not stripped:
                                if opening_comment_lines:
                                    opening_comment_lines.append("")
                                continue
                            if not stripped.startswith("#"):
                                break
                            opening_comment_lines.append(stripped[1:])
                        if "changeable parameters" in " ".join(
                            " ".join(opening_comment_lines).lower().split()
                        ):
                            allowed_blocks.append((block_start, block_end))
                    line_idx = scan_idx + 1
                    break
                scan_idx += 1
            else:
                break

        analyzer = ModuleAnalyzer()
        analyzer.visit(tree)

        for (
            name,
            definition_lines,
        ) in analyzer.module_variable_definitions.items():
            use_lines = analyzer.module_variable_uses.get(name, [])
            if len(definition_lines) != 1 or len(use_lines) != 1:
                continue
            definition_line = definition_lines[0]
            # Unpacking is a normal way to pull out the relevant piece from a
            # function that returns more than one value, so don't warn on
            # single-use names introduced by tuple/list unpacking.
            if (name, definition_line) in analyzer.multi_name_lhs_definitions:
                continue
            # Context-manager aliases are often just resource plumbing between
            # multiple with-items or into the body, so allow single-use names
            # introduced by ``with ... as name``.
            if (name, definition_line) in analyzer.with_alias_definitions:
                continue
            # Loop targets are control-flow variables rather than values to
            # inline, so don't warn on names introduced by ``for ... in ...``.
            if (name, definition_line) in analyzer.for_target_definitions:
                continue
            if any(
                start <= definition_line <= end
                for start, end in allowed_blocks
            ):
                continue
            violations.append(
                Violation(
                    path=str(path),
                    line=definition_line,
                    name=name,
                    kind="module variable",
                    definition_line=definition_line,
                    use_line=use_lines[0],
                )
            )

        for (
            name,
            definition_lines,
        ) in analyzer.top_level_function_definitions.items():
            call_lines = analyzer.direct_function_calls.get(name, [])
            if len(definition_lines) != 1 or len(call_lines) != 1:
                continue
            if (
                name == "main"
                and call_lines
                == analyzer.main_guard_direct_calls.get(name, [])
            ):
                continue
            definition_line = definition_lines[0]
            violations.append(
                Violation(
                    path=str(path),
                    line=definition_line,
                    name=name,
                    kind="top-level function",
                    definition_line=definition_line,
                    use_line=call_lines[0],
                )
            )

    if parse_failures:
        print("\n".join(parse_failures), file=sys.stderr)
        return 1
    if not violations:
        return 0

    report_lines = [
        (
            "Reminder: if a module-level value is intentionally used once "
            "because it is meant to stay easy to edit, place it near the top "
            "of the file inside a vim fold block that starts with '# {{{' and "
            "ends with '# }}}'. The opening comment for that block must "
            "include 'changeable parameters'; that wording may wrap onto "
            "following comment lines."
        ),
        "",
    ]
    for violation in sorted(
        violations, key=lambda item: (item.path, item.line, item.name)
    ):
        if violation.kind == "module variable":
            report_lines.append(
                (
                    f"{violation.path}:{violation.line}: "
                    f"`{violation.name}` is a module-level variable that is "
                    "defined once and used once.\n"
                    f"Defined on line {violation.definition_line}; the only "
                    f"use is on line {violation.use_line}.\n"
                    "The easiest solution is to move the expression from "
                    f"{violation.definition_line} and evaluate in-place on "
                    f"{violation.use_line}\n"
                )
            )
            continue
        report_lines.append(
            (
                f"{violation.path}:{violation.line}: "
                f"`{violation.name}` is a top-level function that is defined "
                "once and directly called once.\n"
                "The easiest solution is to simply move the code to where it "
                "is used.\n"
                "**NOTE** I'm a little more worried about this error message "
                "vs. the variable one.  If you want the function to be "
                "accessible outside the module, you don't want to do this, "
                "and we should discuss ← (JMF)."
            )
        )
    report_lines.append("")
    print("\n\n".join(report_lines))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
