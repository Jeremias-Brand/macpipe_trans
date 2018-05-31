# check if the necessary dirs exist and if not creates them
import os

# 
dirs = ["logs", "logs/trinity", "logs/transrate", "trinity", "transrate"]

def save_mkdir( dirs ):
	for d in dirs:
		if not os.path.isdir(d):
			os.mkdir(d)
<<<<<<< HEAD
			print("Creating directory: " + d)

save_mkdir(dirs)

=======
			print("Creating directory: " + d + ".")

save_mkdir(dirs)
>>>>>>> origin/master
