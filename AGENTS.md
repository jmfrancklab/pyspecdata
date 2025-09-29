# Repository Instructions

- Before running tests or making code changes, always install the project in editable mode with Meson by running `pip install -e . --no-build-isolation` from the repository root.
- Ensure all runtime and testing dependencies are installed so that `pytest` can execute without missing modules. If the editable install fails because build tools are absent (e.g., Meson, ninja, numpy), install the missing packages and rerun the command.
- After installing dependencies, run the relevant `pytest` targets to validate changes.
