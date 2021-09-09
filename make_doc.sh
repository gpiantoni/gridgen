rm docs/gridloc -fr
pdoc3 gridloc -o docs --html
python3 gridloc/bin/parameters.py
git add docs

