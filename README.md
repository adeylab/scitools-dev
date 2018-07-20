
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

      FOLDER -> sci_commands   This folder contains perl module ".pm", python ".py", and R ".r"
                               scitools commands.

To add a new command, do the following:

   1) For a perl command, copy the empty_module.pm file and rename it with the command name.
      The command can only have underscores in the name (will be dashes for calling it). Do
      the same for .py and .r equivalents.

   2) Add in the command text - use the same options structure and help message as other
      commands for consistency.

   3) Test the command and once functional, push it to the master branch. If a number of
      commands or modules will be used (particularly for sci_utils files), creating a new
      branch is recommended.
   
For any command type, once it is developed, do the following:

   1) Within the scitools executable there are two subroutines: load_aliases and load_descriptions

   2) In load_aliases, add another line with the command name and include alternative ways
      to call it if there are any.

   3) In load_descriptions, add another line with a brief statement of the command. This is the
      description that gets listed in "scitools list"
