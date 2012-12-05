#!/bin/bash
for i in {1..20}
do
		# Create the directory
#		mkdir "/Users/apk/Desktop/desktop/startup_weekend_sensor_data/${i}"
    # Check to see if the volume is mounted
    drive="Q${i}"
		if mount|grep $drive;
		then
	    echo "${drive} is mounted"
	    # move the files over to the directory
	    cp -r /Volumes/${drive}/data/12_02_2012 /Users/apk/Dropbox/startup_weekend/analyses/data_working/startup_weekend_sensor_data/${i}/
		else
	    echo "${drive} is NOT mounted"
		fi
done

