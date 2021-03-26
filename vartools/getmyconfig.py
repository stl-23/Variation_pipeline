import configparser
def getConfig(section, key):
    config = configparser.ConfigParser()
    path = os.path.abspath(os.path.dirname(__file__))+'/softwares.config'
    config.read(path)
    return config.get(section, key)