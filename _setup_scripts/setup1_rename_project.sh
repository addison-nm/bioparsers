#!/usr/bin/env bash
#=============================================================================
#
# FILE: rename_project.sh
#
# USAGE: rename_project.sh [name]
#
# DESCRIPTION: Replaces all occurences of the default name "myproject" with the
#   given new name, or prompts for new name. Only handles occurences in the 
#   default project state, and should be run immediately after cloning.
#
# EXAMPLE: sh _setup_scripts/rename_project.sh mynewproject
#=============================================================================

if [ "$#" -gt 1 ]; then
    echo "Usage: $0 [name]"
    exit 1
fi

if [ "$#" -eq 1 ]; then
    newname=$1
else
    read -p "Enter new project name: " newname
fi

# Validate that the name is a valid Python package name
if ! echo "$newname" | grep -qE '^[a-z][a-z0-9_]*$'; then
    echo "Error: '$newname' is not a valid Python package name."
    echo "Name must start with a lowercase letter and contain only lowercase letters, digits, and underscores."
    exit 1
fi

OLDNAME="myproject"

read -p "This action will alter current files and directories. 
Occurences of '$OLDNAME' will be replaced with '$newname'. 
Do you want to continue? (y/n): " answer

case "$answer" in
    [Yy]* ) echo "Continuing...";;
    [Nn]* ) echo "Aborting."; exit 1;;
    * ) echo "Please answer [y]es or [n]o."; exit 1;;
esac


# Replace all occurences in specific files
flist=(
    "pyproject.toml"
    "src/${OLDNAME}/__main__.py"
    "tests/test_core.py"
    "README.md"
)

for f in "${flist[@]}"; do
    echo "Modifying $f"
    sed -i.bak "s/${OLDNAME}/${newname}/g" "${f}" && rm "${f}.bak"
done

# Rename src directory
mv "src/${OLDNAME}" "src/${newname}"
