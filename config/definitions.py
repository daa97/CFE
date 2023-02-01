import os
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))



#Now, use this code whenever you want to navigate to a certain file:
'''import os
from config.definitions import ROOT_DIR
print(os.path.join(ROOT_DIR, 'data', 'mydata.json'))'''