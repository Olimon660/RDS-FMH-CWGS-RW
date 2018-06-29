rm myDatabase
sqlite3 myDatabase
.read ../OTHERS/ANNOTATION.sql
.separator ';'
.import  A.csv  Annotation
