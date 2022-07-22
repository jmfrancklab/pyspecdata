Finding Files
=============

It's a frequent problem that you want to be able to access files on different
instruments, where the absolute path to the file might be in a different
location.
To get around this we define the `exp_type` (experiment type), which points to
the location where a particular type of experiment is stored.
The `exp_type` should have the following characteristics:

-   It should be unique (in a case-insensitive way)
-   It should be the name of the last one to three (separated by a *forward* slash) directories/folders where you have a particular type of experiment stored
-   The name of the folders referred to in `exp_type` should match on the server (google drive, ftp, etc) where the files are stored and on your local computer.

We provide the command XXXX.  You supply this with a csv file containing two columns:

-   `exp_type`
-   file name

Of the files that you want.  It will then run rclone to copy all those files to your local computer.
