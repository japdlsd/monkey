import main

if __name__ == "__main__":
    config = main.defaultConfig()
    with open('defaultConfig.ini', 'w') as f:
        config.write(f)    

