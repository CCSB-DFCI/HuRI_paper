import os
import configparser
import MySQLdb as mysql


def get_connection():
    """Connect to CCSB database.

    If run for first time, asks for username and password and saves
    them in a text file.

    Returns:
        mysql.connect: database connection.
    """
    config = configparser.ConfigParser()
    configPath = '../db_config.ini'
    host = 'paros.dfci.harvard.edu'
    if os.path.exists(configPath):
        config.read(configPath)
        un = config['db']['username']
        pw = config['db']['password']
    else:
        config.add_section('db')
        un = input('Enter your username for database.')
        pw = input('Enter your password for database.')
        config.set('db', 'username', un)
        config.set('db', 'password', pw)
        print('Saving username and password so you dont have to enter them again.')
        with open(configPath, 'w') as f:
            config.write(f)
    return mysql.connect(host, user=un, passwd=pw)


def main():
    get_connection()


if __name__ == '__main__':
    main()
