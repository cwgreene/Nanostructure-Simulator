import sqlite3

parameter ={"parameter_name":"string",
	    "parameter_value":"double",
	    "parameter_units":"string"}
 

tables = {
	#runs
	"runs": 			{"run_id":"id",
		 			 "date":"date",
					 "tag":"string"},
	"runparameters": 		{"run_id":"id",
					 "_param":parameter},
	#frames
 	"framedata": 			{"run_id":"id",
		      			 "timestep":"int",
		      			 "_param":parameter},
	"framepointdata":		{"run_id":"id",
		          		 "timestep":"int",
			  		 "mesh_point_id":"int",
			  		 "_param":parameter},
	#mesh
	"meshes":			{"mesh_id":"id",
		  			 "mesh_name":"string",
		  			 "dimension":"int"},
	"meshruns":			{"run_id":"id",
		    			 "mesh_id":"id"},
	"meshpoints":			{"mesh_id":"id",
		      			 "point_id":"int",
		      			 "x":"double",
		      			 "y":"double",
		      			 "z":"double"},
	#material
	"materials":			{"material_id":"id",
		     			 "material_name":"string"},
	"materialparameters":		{"material_id":"id",
		              		 "_param":parameter},
	#checksums
	"checksums":			{"run_id":"id",
					 "filename":"string",
					 "filehash":"id"}
	}
