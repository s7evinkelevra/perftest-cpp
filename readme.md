## Introduction
This is a short overview of how to setup and run the simulation model from [this thesis](https://drive.google.com/file/d/1XXQ7DSzafJVMXgs--7MrWBanpVhkn5xt/view?usp=sharing) locally and on the HUMMEL HPC cluster of the University of Hamburg. There is a companion repository for the analysis of the output data of this model that can be found [here](https://github.com/s7evinkelevra/sim-analysis). Disclaimer: there might be better ways to do some or all of this, but this is a relatively painless way i found. Building and compiling c++ code has been finicky in my experience, so if something does not work please feel free to reach out. On top of that, i am no expert in any of this and things might not work as i think they do.

The simulation model is written in c++ (version c++17) with cmake as build tool and gcc as compiler. All the source code is available from a public github repository [here](https://github.com/s7evinkelevra/perftest-cpp). C++ is a compiled language, thus needing a compiler to compile the source code (.h and .cpp files in the project) into an executable file. The executable is specific to the system it is compiled on, so you need to compile it yourself to run it on your system. Cmake is a build system i used, that aims to simplify this process by generating the correct build files for your system. The compiler i used is gcc and is available for unix-like systems (linux on the HUMMEL cluster or macos). Compiling and running the simulation on windows can be archived by using a virtual machine or (as i did, and would recommend) using [WSL](https://learn.microsoft.com/en-us/windows/wsl/install). There are also a host of different ways to compile and/or run the model on windows, but i did not explore them.

For all platforms, the process looks roughly like this:
1. Download the source code.
2. Run cmake to generate the build files for your system.
3. Compile the source code using those build files with the compiler.
4. Run the compiled executable.

## Project Structure
Both the .idea and .vscode folder contain config files for code editors (jetbrains clion and vscode) and are not relevant for this.
The CMakeList.txt file is the configuration file for cmake, .gitignore contains all the files that are not be uploaded to the public repository.
Main.cpp is the main file and the entry into the simulation code, which is contained in the src folder. If you want to start to dig into the code, main.cpp is the place to start.
test.cpp only contains some small performance-testing code and is not relevant for the simulation itself.

## Prerequisites
### Windows
- Install WSL and your linux distribution of choice (i recommend ubuntu). There are many resources available on how to enable it, i used [this](https://learn.microsoft.com/en-us/windows/wsl/install). [A virtual machine](https://learn.microsoft.com/en-us/virtualization/hyper-v-on-windows/quick-start/quick-create-virtual-machine) also works.
- Do everything mentioned here in WSL or the virtual machine.

### MacOS arm cpu (M1/M2 Macbook etc.)
To use openMP on arm cpus, you can install gcc through [homebrew](https://brew.sh/) (package manager for macos). This uses the x86 version of gcc through Rosetta2, which supports openMP. You can also install cmake via homebrew. When generating the build files, it is important to specify the homebrew version of gcc to use for compiling (see [[How to Sim#Generating the build files#Macos arm cpu:|this section]]). After that, everything works pretty much the same way.

### MacOS x86 cpu
Works the same as on linux.

### All Platforms
- Local copy of the code
  - If you have git installed, you can use `git clone https://github.com/s7evinkelevra/perftest-cpp` to download the entire repository
  - You can also download it directly from the github page as zip
  - On Hummel, i recommend to use the git module to clone the repository. But you can also upload the files directly using SSH or FTP
- cmake, minimum version 3.15
- gcc, i used 10.3.0 but lower versions also appear to work

If you're running locally, please make sure to have both cmake and gcc installed. If you're running on HUMMEL, you need to load these as "modules". This can be done using `module load [name]`
To load the correct modules, use these commands:

```bash
module load git
module unload env/system-gcc
module load env/gcc-10.3.0
module load cmake
```

The `module unload env/system-gcc` is used to unload the default gcc compiler before loading a newer version. I recommend adding these commands to the .profile file in your user folder to load them automatically when you connect to hummel.

## Generating the build files
Once the the source code is downloaded and everything is installed on the system, we can generate the build files using cmake.

### All Platforms
1. Navigate to the root folder of the downloaded source code.
2. Run cmake using the following commands:
   `cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release`
   This will run cmake with the current directory as the source code (`-S .`), the `./build` folder as output folder for the build files (`-B ./build`) and use the `Release` mode (`-DCMAKE_BUILD_TYPE=Release`). In the Release mode, all debugging information is removed and high optimizations are to be used when compiling.
   If this command ran successfully, you should find you build files in the `/build` folder in the project directory.

### MacOS arm cpu
As mentioned, you need to specify the x86 version of the compiler when generating the build files. This can be done by setting two environment variables when running cmake:
`CC=<path to homebrew gcc> CXX=<path to homebrew g++> cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release`
You can see the path of you homebrew install path by running `brew ls gcc`.  The full command should then look something like this:
`CC=/opt/homebrew/Cellar/gcc/12.2.0/bin/gcc-12 CXX=/opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12 cmake -S . -B ./build DCMAKE_BUILD_TYPE=Release`

## Compiling
To compile the source code using the build files, navigate to the `/build` folder and run
`cmake --build .`
Depending on your system, this might take a few moments. After compiling, you should see the executable `perftest_cpp` in the folder.

## Running the model
When running the model, it looks for a configuration file (`config.json`) either at the path you supplied when running the model as parameter or in the same folder as the executable. This config specifies all the parameters used in the model. For a full explaination of the parameters, please take a look in the methods part of the [thesis](https://drive.google.com/file/d/1XXQ7DSzafJVMXgs--7MrWBanpVhkn5xt/view?usp=sharing). The build tool is setup to copy an example version to the build folder, so there should be one present already. To run the model, simply type `./perftest_cpp`
This will run the model with the default config, producing a number of different output files, located in the `/output` folder. If you want to provide a config file at a different location, run `./perftest_cpp <path-to-your-config>`.

### On Hummel
To utilize the capacities of hummel, we can use [SLURM](https://slurm.schedmd.com/documentation.html) to schedule and run many instances of the simulation in parallel.
__IMPORTANT:__ on hummel, processes ran from a slurm job can only read from the `$HOME` directory, you need to be in the `$WORK` directory to write files. So when executing any of the slurm commands, make sure you `cd $WORK/<your-work-dir>` before.

You can define SLURM jobs using a bash script, with many different options to setup your job using the simulation model. Here, i will briefly go through the script i used to feed each instance a different configuration file. It requires a certain setup, but worked very nicely for my purposes.

For it to work, you need the following folders and files:
- The script below somehere in your `$HOME` folder.
- The executable simulation model at `$HOME/perftest-cpp/build/perftest_cpp`
- A `$HOME/config-queue/` folder that contains the configs you want to run the model with. The config files should be numbered 1 to [number of instances you want to run] + the json file extension. For 6 instances (like below), the config folder thus should contain 6 config files named `1.json, 2.json, 3.json...6.json`

To run it, need to go to your `$WORK` dir and then run `$HOME/<script_name>`.

```bash
#!/bin/bash

######
# run from $WORK -> slurm tasks can only write files there
######

# Set the name of the job
#SBATCH --job-name=simjob
# Specify the partition to use. See https://www.rrz.uni-hamburg.de/services/hpc/hummel-2015/batch-verarbeitung/partitionen-und-job-limits.html (german only...)
#SBATCH --partition=std
# The number of tasks to start, this is equivalent to the number of individual instances started by slurm
#SBATCH --ntasks=6
# How many nodes to occupy. A node is the "computation unit", with one or multiple cpus, ram, etc.
#SBATCH --nodes=6
# How many cpu cores per task/individual instance of the model to use. A single instance can use many cpu cores to speed up the computation significantly. 
# A single std partition node has 2 cpus with 8 cores each (see https://www.rrz.uni-hamburg.de/services/hpc/hummel-2015/hardware-konfiguration.html#knotentypen).
# You may notice that we request 32 cpus per task/node, while only 16 are available. For some reason, it helps with cpu utilization to run the model on 32 threads while only 16 cores are avaible. I'm not entirely sure why, maybe has to do with hyperthreading.
#SBATCH --cpus-per-task=32
# The maximum runtime of a single task. For regular users, it is capped to 12 hours maximum, but you can request longer runtimes or request access to the stl partition (partition for long running jobs).
#SBATCH --time=12:00:00
# Can't remember what this does
#SBATCH --export=NONE

# initialization
source /sw/batch/init.sh

# Set the number of threads used by the simulation model, should be equal to the number f --cpus-per-task
export OMP_NUM_THREADS=32                # essential
#export OMP_PROC_BIND=true               # recommended
#export OMP_PLACES=cores                 # recommended
#export OMP_SCHEDULE=static              # recommended
#export OMP_DISPLAY_ENV=verbose          # good to know
#export KMP_AFFINITY=verbose             # good to know

# $SLURM_PROCID goes from 1 to number of tasks
# e.g. for task 1, it is set to 1, for task 4 to 4, etc.
echo "process id: $SLURM_PROCID"
echo "pwd: $(pwd)"


echo "copying run configs files from config queue"
# remove old configs
rm -r ./configs/.
# copy run configs to working directory
cp -r $HOME/config-queue/. ./configs/

# remove old binary
rm ./perftest_cpp
# copy the current binary to work dir
cp $HOME/perftest-cpp/build/perftest_cpp ./
# make binary executalbe
chmod +x ./perftest_cpp

# run the sim model with srun
# slurm will insert the appropriate $SLURM_PROCID when running the command
# e.g. task number 1 will run "./perftest_cpp ./configs/1.json" and task 4 will run "./perftest_cpp ./configs/4.json" etc.
srun bash -c './perftest_cpp ./configs/$SLURM_PROCID.json'
```

Each instance of the model will produce a folder in `./output` with the config id as name with the simulation output files (the config id is specified in the config.json file).