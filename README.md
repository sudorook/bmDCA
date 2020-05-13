# Boltzmann-machine Direct Coupling Analysis (bmDCA)

Dependencies (installation instructions detailed below):
 * [GCC](https://gcc.gnu.org/) that supports the C++11 standard and
   [OpenMP](https://en.wikipedia.org/wiki/OpenMP)
 * AutoTools
 * pkg-config
 * [Armadillo](http://arma.sourceforge.net/)

This repository contains a C++ reimplementation of bmDCA adapted from [the
original](https://github.com/matteofigliuzzi/bmDCA) code. Method is described
in:

> Figliuzzi, M., Barrat-Charlaix, P. & Weigt, M. How Pairwise Coevolutionary
> Models Capture the Collective Residue Variability in Proteins? Molecular
> Biology and Evolution 35, 1018–1027 (2018).

This code is designed to eliminate the original's excessive file I/O and to
parallelize the MCMC in the inference loop.

## Installing dependencies

__GCC__ is used to compile the source code (and dependencies, if necessary).
The code relies on the `fopenmp` flag for parallelization, so GCC is preferred
over Clang. It also needs support for the C++11 standard, so any GCC later than
version 4.2 will suffice.

__AutoTools__ are a set of programs used to generate makefiles for
cross-platform compilation and installation.

__pkg-config__ is a program that provides a simple interface between installed
programs (e.g. libraries and header files) and the compiler. It's used by
AutoTools to check for dependencies before compilation.

__Armadillo__ is a C++ linear algebra library. It's used for storing data in
matrix structures and performing quick computations in the bmDCA inference
loop. To install, again look to your package repository.

### Linux

To install the dependencies in Linux, simply use your distributions package
manager. Commands for Debian/Ubuntu and Arch Linux are provided below:

#### Debian/Ubuntu

Run:
```
sudo apt-get update
sudo apt-get install git gcc g++ automake autoconf pkg-config \
  libarmadillo-dev libopenblas-dev libarpack++2-dev
```

#### Arch Linux

For Arch Linux, GCC should have been installed with the `base` and `base-devel`
metapackages (`sudo pacman -S base base-devel`), but if not installed, run:
```
sudo pacman -S gcc automake autoconf pkgconf
```

For Arch, Armadillo is not in the package repositories. You will need to check
the AUR.

First, install the SuperLU library:
```
git clone https://aur.archlinux.org/superlu.git
cd superlu
makepkg -si
cd ..
```

SuperLU is a fast matrix factorization library required as a build dependency
for Armadillo. Other build dependencies will be installed via `makepkg` from
the official repositories.

Now, download and install Armadillo:
```
git clone https://aur.archlinux.org/armadillo.git
cd armadillo
makepkg -si
cd ..
```

<!-- If there is no package for Armadillo, or you do not have root privileges on the
   - system your using, you can instead compile the library from source.
   - 
   - First, make sure that `cmake`, `openblas` (or `blas`), `lapack`, `arpack`, and
   - `SuperLU` are installed. CMake is a compilation tool and the others are build
   - dependencies. Then, to download and install Armadillo system wide, run the
   - following:
   - ```
   - wget https://sourceforge.net/projects/arma/files/armadillo-9.850.1.tar.xz
   - tar xf armadillo-9.850.1.tar.xz
   - cd armadillo-9.850.1
   - cmake .
   - make -j4
   - sudo make install
   - cd ..
   - ``` -->


### macOS

The macOS instructions rely on Xcode developer tools and Homebrew for package
management. All commands will be entered into the Terminal.

First, install Xcode developer tools. Open the 'Terminal' application from the
launcher and run:
```
xcode-select --install
```

This may already have this installed.

Next, install Homebrew. From the [online instructions](https://brew.sh), run:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

If you run into permissions errors when installing Homebrew, complaining that
root owns the `/usr/local/` directory, you can change the ownership by running:
```
sudo chown -R <user> /usr/local/
```

where `<user>` should be substituted with your username, e.g. `john`.

Once Homebrew is installed, run:
```
brew install gcc automake autoconf pkg-config armadillo
```

This will install the most recent GCC (9.3.0 as of writing) along with
AutoTools and pkg-config.

__IMPORTANT:__ The default `gcc`, located in `/usr/bin/gcc` is actually aliased
to `clang`, which is another compiler. While in principle this is not an issue,
this version of Clang is not compatible with the `fopenmp` compiler flag that
is used to enable parallelization of the MCMC sampler. Additionally, libraries
(see Armadillo in the next step) installed via Homebrew are not by default
known to `pkg-config` or the linker.

Addressing all of these issues involves overriding the `CC` and `CXX`
environmental variables with the new GCC, updating `PKG_CONFIG_PATH` with paths
to any relevant \*.pc files, and updating `LD_LIBRARY_PATH` with any shared
object library linked at compile time.

Doing this for the first time is a bit bewildering, so for convenience, use the
`rcparams` file in the `tools` directory in this repository. In it are a few
helper functions and aliases. Simply append the contents of that file to your
shell run commands. If you don't know what shell you're using, run:
```
echo $SHELL
```

For bash, copy the contents of `rcparams` to `${HOME}/.bashrc`, and for zsh,
copy to `${HOME}/.zshrc`. The general idea is that macOS versions <=10.14
(Mojave and earlier), uses bash as the default shell, and for >=10.15 (Catalina
and later), Apple switched the default shell to zsh.

You can append the `rcparams` file by copy-pasting the code in your favorite
text editor. You could also do something like `cat tools/rcparams >>
${HOME}/.bashrc`, for example.

As a side note, your macOS may not actually source the `.bashrc` file by
default. If you notice that adding the `rcparams` function has not effect in
new terminals, check that the `${HOME}/.bash_profile` file exists. In it, there
should be a line like `[ -f $HOME/.bashrc ] && . $HOME/.bashrc`. (If the
`.bashrc` file exists, use `source` on it.) If no such like is there, add it
and reload your terminal.

The libraries and headers will be found via the `pkgconfig_find()` and
`ld_lib_add()` functions specified in the `rcparams` file.

__Note:__ Run commands are executed when the shell starts, not when the files
are edited. To update your shell to reflect changes, you can either run:
```
source ${HOME}/.bashrc
```

Or simply open a new shell. (For remote systems, you can just log out and log
in again.)

<!-- The files will be installed to `/usr/local/include` and `/usr/local/lib` by
   - default. This requires root privileges (hence, the `sudo make install` at the
   - end). If you want to install elsewhere, adjust the above commands:
   - ```
   - wget https://sourceforge.net/projects/arma/files/armadillo-9.850.1.tar.xz
   - tar xf armadillo-9.850.1.tar.xz
   - cd armadillo-9.850.1
   - cmake . -DCMAKE_INSTALL_PREFIX:PATH=<alternate_path>
   - make -j4
   - make install
   - cd ..
   - ```
   - 
   - Here, change `<alternate_path>` to wherever you want, for example `${HOME}` or
   - `${HOME}/.local`. -->

### Windows

Before starting, install [MSYS2](https://www.msys2.org). This program is a
package distribution for GNU/Unix tools that can be used to build programs for
Windows.

The installer defaults work fine, and if prompted, open the "MSYS2" shell in
the dialog window.

Once MSYS2 is installed and open, update the base libraries by running:
```
pacman -Syu
```

This will download and install some packages. You will then be prompted to
close the terminal. Close it and open it again. Then, again run:
```
pacman -Syu
```

This will upgrade the remaining packages packaged in the installer to their
most recent versions.

Next, install the dependencies for bmDCA:
```
pacman -S nano vim git \
   autoconf automake-wrapper pkg-config make \
   mingw-w64-x86_64-toolchain \
   mingw-w64-x86_64-openmp \
   mingw-w64-x86_64-arpack \
   mingw-w64-x86_64-lapack \
   mingw-w64-x86_64-openblas \
   mingw-w64-x86_64-armadillo
```

The above command will installed the required programs in the `/mingw64/bin`
directory. Unfortunately, this directory is not on the default PATH. You will
need to add it manually.

Open your `.bashrc` file in a text editor (e.g. `vim ~/.bashrc`). Nano and Vim
were installed in the above command block.

Once open, add the line (at the end of the file):
```
export PATH="/mingw64/bin:$PATH"
```

Then, close and open the MSYS2 terminal again.

_Optionally, edit the `/etc/pacman.conf` file. Uncomment the line `#Color` and
add the line `ILoveCandy`. Just a cosmetic flourish for `pacman`._

## Installing bmDCA (all platforms)

Now that all the dependencies have been installed, compile and install bmDCA
globally (default: `/usr/local`) by running:
```
git clone https://github.com/ranganathanlab/bmDCA.git
cd bmDCA
./autogen.sh --prefix=/usr/local && \
make -j4 && \
make install
cd ..
```

Depending on your platform, the `make install` command may fail due to
permissions issues. To remedy this you can either run `sudo make install`
instead, or you can specify a different installation directory that does not
require administrator privileges. The latter option is particularly useful when
working on remote system not under your control.

Should you want to specify a local directory, for example `$HOME/.local`, run:
```
./autogen.sh --prefix=${HOME}/.local && \
make -j4 && \
make install
```

You can replace the value to the right of `--prefix=` with any other path.
Note, that you should check that it is on your system PATH.

In the event you with to uninstall `bmDCA`, simply run `sudo make uninstall`
or `make uninstall` as appropriate.

Test the installation by running in the terminal:
```
bmdca
```

If the installation worked correctly, this will print the usage information,
e.g.:
```
bmdca usage:
(e.g. bmdca -i <input MSA> -r -d <directory> -c <config file>)
  -i: input MSA (FASTA format)
  -d: destination directory
  -r: re-weighting flag
  -n: numerical MSA
  -w: sequence weights
  -c: config file
  -h: print usage (i.e. this message)
  -f: force a restart of the inference loop
```

### Additional notes for Windows 10

Though `bmdca` should now be easily invoked from within MSYS2, one can update
the system PATH variable to make the binaries accessible system-wide, such as
from the command prompt or other terminal emulators. To update the PATH:

1. Type 'env' in the start search bar.
2. Click 'Edit the system environment variables'.
3. Click on 'Environment Variables...' toward the bottom of the window that
   opens.
4. Select 'Path' in one of the two selection windows (either 'User variables'
   or 'System variables' is fine)
5. Once 'Path' is highlighted, click 'Edit...'
6. Enter the `/usr/local/bin` as a new PATH entry. You can either:
   - Click 'New' in the new window and enter the path to `/usr/local/bin` in
     the MSYS2 installation folder (default: `C:\msys64\usr\local\bin`).
   - Click the 'Browse...' button and navigate to the `C:\msys64\usr\local\bin`
     directory.
7. When the new entry is added, click 'OK' on all the opened windows to set all
   the changes. You will need to close and re-open terminals for the changes to
   be reflected.

## Usage

### Inference (`bmdca`)

This step will take an input multiple sequence alignment (MSA) and a config
file specifying learning parameters and options and then run an inference loop
to fit values to a Potts model for amino acid frequencies at positions (Potts
fields) and pairs of frequencies at pairs of positions (Potts couplings).

The command line flags are:
 - `-i`: input MSA, FASTA format
 - `-d`: directory where output files are written
 - `-r`: (_optional_) flag to compute re-weighting coefficients for each
         sequence in the alignment, with the goal to not unduly bias inference
         by highly similar sequences arising from the phylogeny (default:
         `false`)
 - `-c`: (_optional_) config file for bmDCA run hyperparameters, such as
         `example/bmdca.conf`
 - `-t`: threshold for computing default sequence weights (default: `0.8`)
 - `-n`: input MSA, numerical format
 - `-w`: file containing sequence weights
 - `-h`: print usage information
 - `-f`: force a restart of inference loop (i.e., start at step 1)

If `-r` is not specified, each sequence will be equally weighted, and if no
config file is supplied, the run will default to hyperparameters hard-coded in
the `initializeParameters()` function defined in `src/run.cpp`. The default
number of iterations for the Boltzmann machine is 2000.

The mapping from amino acids to integers is defined in the following way. Amino
acids are ordered as in the following string "-ACDEFGHIKLMNPQRSTVWY". They are
then mapped to the integer corresponding to their position in the string, minus
one. The gap symbol is mapped to 0, A is mapped to 1, etc...

__Important:__ The MSA processing function does not handle gaps represented by
'.' characters.

#### Example 1: run from FASTA file

To learn a FASTA-formatted multiple sequence alignment (with re-weighting) and
a config file:

```
bmdca -i <input_alignment.fasta> -d <output_directory> -r -c <config_file.conf>
```

#### Example 2: run from weighted numerical alignment

If you already have a numerically-formatted alignment (gaps are 0) and set of
per-sequence weights, run:
```
bmdca -n <numerical_alignment.txt> -w <sequence_weights.txt>
  -d <output_directory> -c <config_file.conf>
```

#### Example 3: restarting runs

Take, for example, this command:
```
bmdca -i <input_alignment.fasta> -d <output_directory> -r -c <config_file.conf>
```

If the run is stopped before the maximum number of steps is reached, simply
invoke the same command again to restart it:
```
bmdca -i <input_alignment.fasta> -d <output_directory> -r -c <config_file.conf>
```

The inference loop will pick up from where it left off.

**IMPORTANT** To guarantee that inferences loops produce the same results
irrespective of whether they were stopped and restarted or ran continuously,
bmDCA will check that the hyperparameters used previously and presently match.
The only fields that may be adjusted before restarting are:

1. `save_parameters`
2. `step_max` - the max number of steps (increase to continue a loop after it
                ends)
3. `error_max` - the convergence threshold (lower it to continue)

The hyperparameters used by any given run are stored in the `bmdca_params.conf`
file. For boolean configuration options, it is not to specify `1` or `0`
instead of `true` and `false`.

### Sampling (`bmdca_sample`)

Use a Monte-Carlo sampler to draw sequences from the model specified by the
learned parameters.

Run:
```
bmdca_sample -p <parameters.txt> -d <output_directory> \
  -o <output_file.txt> -n <number_of_sequences> \
  -r <number_of_indep_sampling_runs> -c <config_file.conf>
```

If instead, you save parameters in binary format, run:
```
bmdca_sample -p <parameters_h.bin> -P <parameters_J.bin>  \
  -d <output_directory> -o <output_file.txt> \
  -n <number_of_sequences> -r <number_of_indep_sampling_runs> \
  -c <config_file.conf>
```

The command line flags are:
 - `-p`: input parameters, text format __or__ fields (h) parameters file,
         binary format
 - `-P`: (_optional_) couplings (J) parameters file, binary format
 - `-d`: directory where output files are written
 - `-c`: (_optional_) config file for bmDCA run hyperparameters, e.g.
         `example/bmdca.conf`
 - `-o`: name of the output file for the sequences
 - `-n`: number of sequences to sample in each independent run (default: 1000)
 - `-r`: number of independent sequencing runs (default: 10)

Note, you only need to specify the `-p` option if the `bmdca` output is stored
in text files. For binary parameters, which are stored in two `_h_%d.bin` and
`_J_%d.bin` files, pass the fields (h) file to `-p` and couplings (J) file to
`-P`.

### Conversion (`arma2ascii`)

If using the `output_binary=true` flag during the inference step, the output
will be stored in an Armadillo-specific binary format. While this allows for
reproducible outputs in stopped-and-restarted inference runs, the format is not
accessible for other programs. You can use the `arma2ascii` tool to convert
binary-stored outputs to ASCII.

To convert parameters:
```
arma2ascii -p <parameters_h.bin> -P <parameters_J.bin>
```

To convert stats files:
```
arma2ascii -s <MC_stat_file.bin>
```

The output fill be stored in a `.txt` file.

## Configuration file options

Inference and sampling runs can be configured using a text file (see
`example/bmdca.conf`). The fields in the file are as follows:

### [bmDCA]

1. `lambda_reg1` - L2 regularization strength for fields, h (default: 0.01)
2. `lambda_reg2` - L2 regularization strength for couplings, J (default: 0.01)
3. `step_max` - maximum number of iterations for Boltzmann learning process
   (default: 2000)
4. `error_max` - error convergence criterion for stopping (default: 1e-05)
5. `save_parameters` - save parameters every `save_parameters` number of steps
   (default: 100)
6. `save_best_steps` - save steps that yield the lowest RMSE of the gradient
   (default: false)
7. `random_seed` - initial seed for the random number generator (default: 1)
8. `use_reparametrization` - use the re-parametrized model for inference
   (default: true)
9. `epsilon_0_h` - initial learning rate for fields (default: 0.01)
10. `epsilon_0_J` - initial learning rate for couplings (default: 0.001)
11. `adapt_up` - multiple by which to increase Potts (J and h) gradient
    (default: 1.5)
12. `adapt_down` - multiple by which to decrease Potts (J and h) gradient
    (default: 0.6)
13. `min_step_h` - minimum learning rate for h (default: 0.001)
14. `max_step_h` - maximum learning rate for h (default: 2.5)
15. `min_step_J` - minimum learning rate for J (default: 1e-05)
16. `max_step_J_N` - maximum learning rate for J, scaled by effective number of
    sequences (default: 2.5)
17. `error_min_update` - threshold for differences in MSA and MCMC frequencies
    above which parameters (J and h) are updated (default: -1)
18. `t_wait_0` - initial burn-in time (default: 10000)
19. `delta_t_0` - initial wait time between sampling sequences (default: 100)
20. `check_ergo` - flag to check MCMC sample energies and autocorrelations,
    without which wait and burn-in times are not updated (default: true)
21. `adapt_up_time` - multiple to increase MCMC wait/burn-in time (default: 1.5)
22. `adapt_down_time` - multiple to decrease MCMC wait/burn-in time (default
    0.6)
23. `step_important_max` - maximum number of importance sampling steps
    (default: 1, i.e.importance sampling disabled)
24. `coherence_min` - (default=.9999)
25. `use_ss` - flag to sample sequences equal to the effective number of
    sequences in the alignment (default: false)
26. `M` - number of sequences to sample for each MCMC replicate (default: 1000)
27. `count_max` - number of independent MCMC replicates (default: 10)
28. `init_sample` - flag for whether of not to use seed sequence for
    initializing the MCMC (default: false)
29. `init_sample_file` - file containing the MCMC seed sequences (default: "")
30. `sampler` - sampler mode, 'mh' for Metropolis-Hastings and 'z-sqrt' or
    'z-barker' for Zanella, 2019. 'z-sqrt' corresponds to a balancing function of
    `sqrt(t)`, and 'z-barker' corresponds to `t/(1+t)`. (default: "mh")
31. `use_pos_reg` - flag to apply position-specific regularization when
    learning J (default: false)
32. `output_binary` - flag to output data in binary format, which is faster and
    more precise (default: true)

### [sampling]

1. `random_seed` - initial seed for the random number generator (default: 1)
2. `t_wait_0` - initial burn-in time (default: 100000)
3. `delta_t_0` - initial wait time between sampling sequences (default: 1000)
4. `check_ergo` - flag to check MCMC sample energies and autocorrelations,
   without which wait and burn-in times are not updated (default: true)
5. `adapt_up_time` - multiple to increase MCMC wait/burn-in time (default: 1.5)
6. `adapt_down_time` - multiple to decrease MCMC wait/burn-in time (default
   0.6)
7. `sampler` - sampler mode, 'mh' for Metropolis-Hastings and 'z-sqrt' or
   'z-barker' for Zanella, 2019. 'z-sqrt' corresponds to a balancing function
   of `sqrt(t)`, and 'z-barker' corresponds to `t/(1+t)`. (default: "mh")
8. `temperature` - temperature at which to sample sequences (default: 1.0)

## Output files

`bmdca` will output files during the course of its run:
 - `bmdca_params.conf`: a list of the hyperparameters used in the learning
   procedure.
 - `energy_%d.dat`: mean and std dev over replicates for sample sequence
   energies at each step of the Markov chain
 - `ergo_%d.dat`: set of autocorrelation calculations for sampled sequences
   used for deciding whether to increase/decrease MCMC wait intervals and
   burn-in times
   1. correlation of sequences 1 wait interval apart.
   2. correlation of sequences M/10 wait intervals apart. (M = # sequences)
   3. cross correlation of sequences
   4. standard deviation of correlations 1 wait intervals apart
   5. standard deviation of correlations M/10 intervals apart
   6. standard deviation of cross correlations
   7. combined deviation of cross and autocorrelations (1 wait interval)
   8. combined deviation of cross and autocorrelations (M/10 wait intervals)
   9. combined deviation of autocorrelations 1 and M/10 intervals apart
 - `MC_energies_%d.txt`: energies of each MCMC sequence, grouped by replicate
 - `MC_samples_%d.txt`: sequences sampled from MCMC, grouped by replicate
 - `msa_numerical.txt`: numerical representation on input MSA
 - `my_energies_cfr_%d.txt`: statistics of energies over replicates, used for
   deciding whether to increase/decrease MCMC wait intervals and burn-in times:
   1. number of replicates
   2. average over replicates of the energies of starting MCMC sequences
   3. standard deviation over replicates of energies of starting MCMC sequences
   4. number of replicates
   5. average over replicates of the energies of ending MCMC sequences
   6. standard deviation over replicates of energies of sending MCMC sequences
 - `my_energies_cfr_err_%d.txt`: additional energies statistics
   1. average over replicates of the energies of starting MCMC sequences
   2. average over replicates of the energies of ending MCMC sequences
   3. combined std dev of energies for staring and ending MCMC sequences
 - `my_energies_end_%d.txt`: energies of ending MCMC sequence for each replicate
 - `my_energies_start_%d.txt`: energies of starting MCMC sequence for each replicate
 - `overlap_%d.txt`: overlap  of pairs of MCMC sequences
   1. number of steps apart (in units of wait time)
   2. mean overlap for all sequences %d steps apart
   3. standard deviation of overlaps for all sequences %d steps apart
 - `parameters_%d.txt`: learned Potts model parameters (J and h)
 - `parameters_h_%d.bin` and `parameters_J_%d.bin`: learned Potts model
   parameters (J and h), stored in arma binary format (see `output_binary=true`
   flag from config file).
 - `rel_ent_grad_align_1p.txt`: relative entropy gradient for each amino acid
   at each position
 - `sequence_weights.txt`: weights for each sequence, either a number between 0
   and 1 based on sequence similarity or 1 if re-weighting was not specified
 - `stat_align_1p`: table of frequencies for each amino acid at each position
   in the MSA
 - `stat_align_2p`: table of frequencies for pairs of amino acids at each pair
   of positions in the MSA (due to symmetry, only the 'upper triangle' of
   positions is stored)
 - `stat_MC_1p_%d.txt`: table of frequencies for each amino acid at each
   position of the set of MCMC-sampled sequences.
 - `stat_MC_1p_sigma_%d.txt`: table of standard deviation of frequencies over
   replicates for each amino acid at each position of the set of MCMC-sampled
   sequences.
 - `stat_MC_2p_%d.txt`: table of frequencies for pairs of amino acids at each
   pair of positions from the set of MCMC-sampled sequences
 - `stat_MC_2p_sigma_%d.txt`: table of standard deviation over replicates of
   frequencies for pairs of amino acids at each pair of positions from the set
   of MCMC-sampled sequences

The final outputs will be stored with a `_final` suffix in the file name before
the file extension. For example, the final learned parameters will be stored in
`parameters_final.txt`. __Use this file or the latest learned parameters for
sampling synthetic sequences.__

Depending how many times you configure `bmdca` to save steps to disk, the total
data generated can be substantial ( > 1 Gb). At present, the only way to
disable writing of a particular log file is to comment out the code in the
`Sim::run()` function defined in `src/run.cpp`.

## Output file formats

### Numerical sequence alignment

This file is a space-delimited file, e.g.:
```
4914 53 21
0 2 10 10 13 16 1 7 6 13 2 1 12 19 17 17 15 19 20 5 18 6 18 18 6 15 2 12 15 5 19 20 6 6 2 7 6 12 9 12 16 5 1 16 4 4 4 2 11 15 18 2 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 6 20 9 6 2 6 6 12 9 12 12 5 3 17 3 6 17 2 17 16 10 2 9
```

The first line is:
1. Number of sequences (M)
2. Number of positions (N)
3. Size of amino acid alphabet (all AAs + 1 for gaps) (Q)

### Learned Potts model parameters

#### Armadillo binary

By default, `bmdca` will save learned parameters in binary format.  These files
(`parameters_h_%d.bin` and `parameters_J_%d.bin`) cannot be directly viewed by
a text editor. To view the contents, convert the files to ASCII by using the
provided `arma2ascii` tool.

See the above usage section for how to use `arma2ascii`. Parameters converted
by the program will match the format for parameters generated when
`output_binary=false` is specified in the config file. See the below section
for details.

#### ASCII

Learned parameters saved in text files are called `parameters_%d.txt`. They
contain the parameters for both J and h, formatted as follows:
```
J [position index i] [position index j] [amino acid index a] [amino acid index b]
.
.
.
h [position index i] [amino acid index a]
.
.
.
```

The position indices go from 0 to N-1 (N = # positions), and the amino acid
indices go from 0 to 20 (21 amino acids total, including gaps). 0 corresponds
to a gap.

### Sequence statistics

The sequence statistics files (e.g. `stat_align_1p.txt` and
`stat_align_2p.txt`) have a different format.

For 1 position (1p) frequencies:
```
[position index] [amino acid frequencies (21)]
.
.
.
```
where `[amino acid frequencies (21)]` is a row of frequencies for each of the
21 positions.

For 2 position (2p) frequencies:
```
[position index i] [position index j] [amino acid frequencies (21x21)]
.
.
.
```
where `[amino acid frequencies (21x21)]` is a row that corresponds to the
frequencies of the `21x21` pairs of amino acids at positions i and j.

## Extra

__For users of shared resources:__ OpenMP will default to the number of
available cores, so if the bmDCA programs are run on a shared resource, say a
cluster, all cores will be engaged, starving other processes of resources or
getting you booted off the system. To prevent this, use the `OMP_NUM_THREADS`
environmental variable.

You can either set it at runtime:
```
OMP_NUM_THREADS=4 bmdca -i ...
```

Or, you can set it globally, for example in your shell rc file.
```
export OMP_NUM_THREADS=4
```

The above examples will limit OpenMP to 4 threads.

__You don't need to worry about this if submitting jobs through a workload
manager__, such as Slurm or Sun Grid Engine. The manager will limit bmDCA to
the number of cores specified, so manipulating `OMP_NUM_THREADS` is not needed.
