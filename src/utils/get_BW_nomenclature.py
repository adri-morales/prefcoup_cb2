import urllib.request
import re
from bs4 import BeautifulSoup



def get_bw(pdb):

    """
    Takes an url from the protein page of GPCRdb and returns a dictionary linking position and bw nomenclature.

    :param url: URL of the protein page in GPCRdb of structure: https://gpcrdb.org/protein/<protein name>/
    :return: Dictionary with the absolute position as key and the ballesteros-weinstein nomenclature as values
    """
    
    # Get protein name from pdb
    url = f"https://gpcrdb.org/structure/{pdb}"
    fp = urllib.request.urlopen(url)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()
    
    protein = re.findall('\.\./protein/(\w+)"', mystr)[0]

    url =  f"https://gpcrdb.org/protein/{protein}"
    # parsing the html file from the page
    fp = urllib.request.urlopen(url)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()

    # replacing newlines with spaces
    mystr = " ".join(mystr.splitlines())

    # searching strings with position and ballesteros informatino
    matches = re.findall("title='(.+?)'", mystr)

    # spliting strings in two or one element lists
    dic = [title.split() for title in matches]

    # making a dictionary out of the list. residues without BWid will have the same string as key and value
    dic = dict([a if len(a) > 1 else a*2 for a in dic])

    return dic


if __name__ == '__main__':
    d = get_bw('7wu2')
    print(d)