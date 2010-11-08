import database_schema as dbs
import sqlite3
import optparse
import sys

types = {"id":"text",
	 "string":"text",
	 "date":"text",
	 "double":"real",
	 "int":"integer"};

def unroll_attr(complex,unrolled):
	result = []
	for attr in complex:
		if attr[0] == "_":
			unroll_attr(complex[attr],unrolled)
		else:
			unrolled[attr]=types[complex[attr]]
	return 

def create_table(name,attributes):
	statement = "create table %(name)s (%(attrlist)s);"
	attrlist = ""
	unrolled={}
	unroll_attr(attributes,unrolled)
	for attr in unrolled:
		attrlist+= attr+" "+unrolled[attr]+","
	attrlist = attrlist[:-1] #get rid of trailing comma
	statement = statement %{"name":name,"attrlist":attrlist}
	return statement

def init(filename):
	database = sqlite3.connect(filename)
	tables = dbs.tables
	statement = ""
	for table in dbs.tables:
		database.execute(create_table(table,tables[table]))

if __name__=="__main__":
	init(sys.argv[1])
