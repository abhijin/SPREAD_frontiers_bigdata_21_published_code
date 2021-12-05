import os
import io
import sqlite3
from sqlite3 import Error
import csv
import matplotlib.pyplot as plt
#import scipy.stats as stats

def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn


def update_task(conn, task):
    """
    update priority, begin_date, and end date of a task
    :param conn:
    :param task:
    :return: project id
    """
    #sql variable is defunct
    sql = ''' UPDATE projects
              SET   localityType = ?,
                    diameter = ?,
                    mooreRange = ?,
                    rowLength = ?,
                    squareRowSize = ?,
                    whiteNodeRows = ?,
                    LDGraphType = ?,
                    maxCoreNumber = ?,
                    spectralRadius = ?'''
    #inserting into which columns in the sql file order of variables and columns.
    sqlite_insert_with_param = """INSERT INTO projects
                          (localityType, diameter, mooreRange, rowLength, squareRowSize, whiteNodeRows, LDGraphType, maxCoreNumber, spectralRadius)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);"""
    #inserting into sql file.
    cur = conn.cursor()
    cur.execute(sqlite_insert_with_param, task)
    conn.commit()

def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        #Creates a table if table is not already created.
        c = conn.cursor()
        c.execute(create_table_sql)
        conn.commit()
        c.close()
    except Error as e:
        print(e)
def delete_eigen_value():
    """Code
    for getting rid of the imaginary component of the eigen values in my data text document.Better Data without imaginary component
    stored in completeData.txt raw data stored in dataStoring.txt
    """
    dataFile = open("dataStoring.txt", "r")
    fileLines = dataFile.readlines()
    dataFile.close()
    newDataFile = open("completeData.txt", "a")
    for lines in fileLines:
        adjustedLine = lines[:len(lines) - 9]
        adjustedLine += "),2);\n"
        newDataFile.write(adjustedLine)

    newDataFile.close()

def creating_database():
    # finding the absolute path to the directory the program is in
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    my_file = os.path.join(THIS_FOLDER, 'simulationData.db')
    # absolute path
    database = my_file

    #create a database connection
    conn = create_connection(database)

    # name of the table and all of the column names in the table
    sql_create_projects_table = """ CREATE TABLE IF NOT EXISTS projects (
                                        number_of_nodes integer,
                                        number_of_regions integer,
                                        range float,
                                        locality_size integer,
                                        locality_graph text,
                                        long_distance_type text,
                                        max_core_number integer,
                                        spectral_radius float,
                                        diameter integer
                                        ); """
    # create a database connection
    """""
    code for getting the data from the text file into an sql database
    conn = create_connection(database)
    """
        # create tables
    if conn is not None:
        # create projects table
        print("Creating Table")
        create_table(conn, sql_create_projects_table)
    else:
        print("Error! cannot create the database connection.")


    #add all the data from the measures and properties of the graph into the SQL file.
    #update_task(conn, (whichLocality, graphDiameter, Mrange, sideLength, SquareSide, whiteRows, whichLD, maxCoreNum, eigenValue[0]))


    #query = "INSERT INTO projects (spectralRadius) VALUES (1);"
    cur = conn.cursor()

    dataFile = open("completeData.txt", "r")
    fileLines = dataFile.readlines()
    dataFile.close()
    #newDataFile = open("completeData.txt", "a")
    for lines in fileLines:
        #print(lines)
        cur.execute(lines)
        conn.commit()



def plot_eigenvalue_against_range():
    """
    This is used to read all of the eigenvalues from the database and then sort them by range and find the average. 
    """
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    my_file = os.path.join(THIS_FOLDER, 'simulationData.db')
    # absolute path
    database = my_file
    # create a database connection
    conn = create_connection(database)
    cur = conn.cursor()

    cur.execute('''SELECT * from projects''')
    rows = cur.fetchall() #selects all the rows
    eigenValues = [0, 0, 0, 0, 0, 0, 0, ]
    timesAppeared = [0, 0, 0, 0, 0, 0, 0, ]
    avergaeEigenValues = [0, 0, 0, 0, 0, 0, 0, ]
    range = [0, 0, 0, 0, 0, 0, 0, ]
    # print(rows)
    for row in rows:
        index = int((row[2] - 1) / 0.5) #index corresponding to rnage
        eigenValues[index] += row[7] #add all eigen values of same range togehter
        timesAppeared[index] += 1 #counts times each range appeared
        range[index] = row[2] #which index cooresponds to which range

    # print(eigenValues)
    x = 0
    while x < 7:
        avergaeEigenValues[x] = eigenValues[x] / timesAppeared[x] #calculates average eigenvalue
        print(avergaeEigenValues[x])
        x += 1

    plt.plot(range, avergaeEigenValues) #plots range and averageEigenValues

    # naming the x axis
    plt.xlabel('Range')
    # naming the y axis
    plt.ylabel('Spectral Radius')

    plt.title('Spectral Radius vs Range')
    plt.show()

    # f, p = stats.f_oneway(avergaeEigenValues, range)
    # print(f, p)


