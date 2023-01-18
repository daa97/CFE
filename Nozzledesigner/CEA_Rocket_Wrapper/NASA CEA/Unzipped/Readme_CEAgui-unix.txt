/** 
    File: Readme_CEAgui-unix.txt 	
    The CEA GUI Download Readme file for Unix (Sun, Sgi, Linux, Mac OS X). 
    Last update:  September 28, 2005, Minna.M.Chao@nasa.gov  (216) 433-5197	
    http://www.grc.nasa.gov/WWW/CEWeb/ceaguiDownload-unix.htm
 */
	
A) The CEA GUI Download consists of the following files:

   1. CEAgui-jar.tar.Z 	- CEAgui JAR file:
      If you have previously installed CEAgui Package , then update to a 
      NEW version of CEAgui by downloading the newest JAR file.

   2. CEA+Fortran.tar.Z	- CEA Fortran Source code, Libraries, and 14 
      examples files. If the Fortran files (cea2.f, cea.inc, syntax.f, b1b2b3.f)
      must be recompiled on your local machine,then the Fortran executable 
      files for cea2.f, syntax.f, b1b2b3.f should be renamed as FCEA2, syntax, 
      b1b2b3 respectively. It also included CEA related packages such Cap, PAC, and MISC
	
   3. CEAgui Executer Packages are platform dependence:
        (For your platform where xyz is either sun, sgi, or linux )
        CEAexec-xyz.tar.Z - The CEA Executer package consists of the following 
                            executable,library files, and shell scripts:
      	b1b2b3 and syntax               - The executables for CEAgui. 
      	FCEA2                                - The cea2 executable.
      	thermo.lib and trans.lib         - cea2 Fortran library files.  
      	Readme_CEAgui-unix.txt	- The file you are reading.
      	ceaguiDownload-unix.htm     - The Web page for download.
      	runCEA.sh                          - run CEAgui using this shell script 

      		     The shell script must contain the following line:
		         java -classpath CEAgui.jar CEAgui

   4. If your local system has NOT installed the Java2 SDK (J2SDK), you must 
      install the Java 2 Runtime Environment(JRE)which consists of Java Virtual 
      Machine, the Java plateform core classes, and supporting files. It is the
      runtime part of the Java2 SDK and can be freely downloaded from Sun 
      Microsoft Inc.
      
      For Examples the download from 
        "http://java.sun.com/j2se/1.4.2/download.html"
        and for other platforms from http://java.sun.com/cgi-bin/java-ports.cgi:     
        (a) for Sun Solaris/Sparc:
      		j2re-1_4_1_02-solaris-sparc.sh 
        (b) for Linux:
      		j2re-1_4_1_02-linux-i586.bin (self-extracting binary file)
        (c) for Sgi Irix:
      		j2re-1_4_1_02-sgi.sh 
        (D) for Mac OS X 10.3.5 and Apple Xcode 1.2 (Development Environment)
		http://www.versiontracker.com/ searching for "Mac OS X"
   
    
B) The CEA GUI Executer installaton procedures for your platform (say, xyz):
						      
   1. Create the Installation Directory on your local machine. (mkdir CEAexec)

   2. Download all three(3)TAR files and the executable files and save into 
      the SAME directory.
		
   3. Using uncompress and tar to extract the Fortran Source Code 
      (CEA+Fortran.tar.Z).
		        zcat CEA+Fortran.tar.Z  | tar xvf -
		
   4. Using uncompress and tar to extract the CEA GUI Executer. 
      For your platform where xyz is either sun, sgi, or linux 
      (CEAexec-xyz.tar.Z).
		        zcat CEAexec-xyz.tar.Z | tar xvf -
		
   5. Using uncompress and tar to extract the CEAgui JAR file 
      (CEAgui-jar.tar.Z).
 		        zcat CEAgui-jar.tar.Z  | tar xvf -
 		
   6. If your selected platform is Sun, then execute the J2RE file 
      (say, j2re-1_4_1_02-solaris-sparc.sh) to install Java Runtime Environment 
	   by using the selected directory (e.g. CEAexec\JRE\)
 
   7. Update Unix environment variables such as JDK_HOME or JRE_HOME, 
      CLASSPATH, and PATH: 
      In csh,the Unix environment variables are modified with the setenv command
 		    setenv CLASSPATH path1:path2
   		    setenv JDK_HOME /usr/jdk1.4
		    setenv PATH $JDK_HOME/bin:${PATH}:.
				   Note That: The DOT . is necessary.

      In sh, the Unix environment variables are modified with these commands:
		    CLASSPATH = path1:path2
		    export CLASSPATH

		    JDK_HOME = /usr/jdk1.4
		    export JDK_HOME

		    PATH = $JDK_HOME/bin:${PATH}:.
			      Note That: The DOT . is necessary.
		    export PATH

      

   8. Setting file permissions for the following executable files	   	        	
		    chmod a+x b1b2b3
		    chmod a+x syntax
		    chmod a+x FCEA2
		    chmod a+x runCEA.sh
   
   9. Execute the shell script (runCEA.sh) to start CEAgui.
      Try examples by selecting "Open Examples" from the "File" Menu and 
      choosing an example. Then select "Run CEA2 Executable" from the 
      "Activity" Menu.	
