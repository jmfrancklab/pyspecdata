# a python replacement for sed, adapted from here: https://stackoverflow.com/questions/12420779/simplest-way-to-get-the-equivalent-of-find-in-python
def listfiles(folder):
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            if filename.lower().endswith('.py'):
                yield os.path.join(root, filename)
for filename in listfiles('~'):
    # rather than printing, do a search based on word boundaries
    # -- '\b' searches for word boundaries
    # use a spreadsheet
    # this is also a good time to develop ability to search in tex for
    # \begin{python} \end{python} blocks (look for notes on TexSoup, etc),
    # not in openproject
    print(filename)
