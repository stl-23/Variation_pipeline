import configparser
def getConfig(section, key):
    config = configparser.ConfigParser()
    path = os.path.abspath('./softwares.config')
    config.read(path)
    return config.get(section, key)