git diff --name-only -z --diff-filter=ACMR devel -- '*.py' | xargs -0 -r -n1 python .github/scripts/check_single_use.py
git diff --name-only -z --diff-filter=ACMR devel -- '*.py' | xargs -0 -r -n1 ~/myflake.sh
