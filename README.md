
The scitools-dev repository is for internal Adey Lab or invited collaborator usage only.

The structure of scitools is the following:

   FOLDER -> PATH/scitools-dev/
  
      FILE -> scitools         This is the main executable, it loads modules and searches for
                               available commands / scripts and executes them

      FILE -> scitools.cfg     This is the config file for scitools, it is where default
                               files, locations, and external program calls are present.
                               To use your own scitools.cfg, copy it to your home directory.

      FILE -> sci_modes.cfg    This is a read mode specification file to be used with the newest
                               version of fastq-dump.

      FILE -> README.md        This file.

      FILE -> SCI_Indexes.txt  Soon to be considered deprocated, old index file.

      FOLDER -> index_files    This folder contains index files to be used with read_indexdir.

      FOLDER -> sci_utils      This folder contains files / perl modules to be loaded by
                               scitools functions, most notably: general.pm

      FOLDER -> sci_commands   This folder contains perl module ".pm" scitools commands.

      FOLDER -> sci_python     This folder contains python scripts ".py" to be called.

To add a new perl command, do the following:

   1) For a perl command, copy the empty_module.pm file and rename it with the command name.
      The command can only have underscores in the name (will be dashes for calling it).

   2) Add in the command text - use the same GetOpts structure and help message as other
      commands for consistency.

   3) Test the command and once functional, push it to the master branch. If a number of
      commands or modules will be used (particularly for sci_utils files), creating a new
      branch is recommended.

To add a new python script, do the following:

   1) Create the python script in the sci_python folder, the same naming rules apply, as this
      is how the command directory and calling is carried out.

   2) Conform to the same options and help text format as with other perl modules.

   3) scitools will call the command and pass it all arguments (other than the command name)
