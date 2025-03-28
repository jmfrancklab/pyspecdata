#!/bin/bash
FAIL=0  # Initialize the fail variable

for file in $CHANGED_FILES; do
    if grep -q "TODO ☐" "$file"; then
        echo "❌ Found TODO ☐ in $file"
        FAIL=1  # Set fail flag if a TODO is found
    fi
done

if [ $FAIL -eq 1 ]; then
    echo "❌ TODO check failed."
    exit 1
else
    echo "✅ No TODO ☐ found in changed files."
    exit 0
fi
