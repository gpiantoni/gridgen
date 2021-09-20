rm docs/gridgen -fr
pdoc3 gridgen -o docs --html
python3 gridgen/bin/parameters.py
git add docs

