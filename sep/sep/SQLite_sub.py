import sqlite3

def get_prop_of_name(_name):
    try:
        sqlite_connection = sqlite3.connect('db.db')
        cursor = sqlite_connection.cursor()
        #print("База данных создана и успешно подключена к SQLite")

        sqlite_select_query = "select * from test where TaggedName = '" + _name + "'"
        cursor.execute(sqlite_select_query)
        record = list(cursor.fetchall())
        res = record[0]
        sqlite_select_query = "pragma table_info(test)"
        cursor.execute(sqlite_select_query)
        record1 = list(cursor.fetchall())
        res1 = record1
        res_dic = {res1[i][1]: res[i] for i in range(len(res))}

        cursor.close()
        return res_dic


    except sqlite3.Error as error:
        print("Ошибка при подключении к sqlite", error)
    finally:
        if (sqlite_connection):
            sqlite_connection.close()
