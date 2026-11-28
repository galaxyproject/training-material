```
(.venv) rivendell$ sqlite3 database/universe.sqlite 
SQLite version 3.34.1 2021-01-20 14:10:07
Enter ".help" for usage hints.
sqlite> .schema user_favorite_extension 
CREATE TABLE user_favorite_extension (
	id INTEGER NOT NULL, 
	user_id INTEGER, 
	value VARCHAR, 
	PRIMARY KEY (id), 
	FOREIGN KEY(user_id) REFERENCES galaxy_user (id)
);
sqlite> select * from alembic_version;
fe1ce5cacb20
d4a650f47a3c
```
