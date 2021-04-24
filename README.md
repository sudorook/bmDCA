# Boltzmann-machine Direct Coupling Analysis (bmDCA)

Dependencies (installation instructions detailed below):
 * [GCC](https://gcc.gnu.org/) that supports the C++11 standard and
   [OpenMP](https://en.wikipedia.org/wiki/OpenMP)
 * AutoTools
 * pkg-config
 * [Armadillo](http://arma.sourceforge.net/)

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
sudo pacman -S gcc automake autoconf pkgconf superlu
```

For Arch, Armadillo is not in the package repositories. You will need to check
the AUR.
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

This may already be installed.

Next, install Homebrew. From the [online instructions](https://brew.sh), run:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

If you run into permissions errors when installing Homebrew, complaining that
root owns the `/usr/local/` directory, you can change the ownership by running:
```
sudo chown -R <user> /usr/local/
```

where `<user>` should be substituted with your username, e.g. `john`. Your
username should be apparent form the command prompt, but if unsure, run
`whoami`.

Once Homebrew is installed, run:
```
brew install gcc automake autoconf pkg-config armadillo
```

This will install the most recent GCC (10.2.0 as of writing) along with
AutoTools and pkg-config.

__IMPORTANT:__ The default `gcc`, located in `/usr/bin/gcc` is actually aliased
to `clang`, which is another compiler. While in principle it can be used to
compile bmDCA, this version of Clang is not compatible with the `fopenmp`
compiler flag that is used to enable parallelization. Additionally, libraries
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

You can append the `rcparams` file by copy-pasting the code in it to run rc
file using your favorite text editor. Another method for appending is to run
`cat tools/rcparams >> ${HOME}/.bashrc` or `cat tools/rcparams >>
${HOME}/.bashrc`, as appropriate.

The `rcparams` file contains a few helper functions, `pkgconfig_find()` and
`ld_lib_add()`, that automate the process of finding where build dependencies
are installed on your system.

__Note:__ Run commands are only executed when the shell starts, _not_ when the
their files are edited. To update your currently-running shell to reflect new
changes, you can either run in the command prompt:
```
source ${HOME}/.bashrc
```

or simply open a new terminal. (For remote systems, you can just log out and log
in again.)

__Bash users:__ Your macOS installation may not actually source the `.bashrc`
file by default. Check that the functions are actually being loaded by running:
```
LC_ALL=C type pkgconfig_find
```

If the function is not found, check that the `${HOME}/.bash_profile` file
exists. In it, there should be a line that looks like `[ -f $HOME/.bashrc ] &&
. $HOME/.bashrc`. If no such like is there or if no such file exists, add it
and open a new terminal.


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

This will download and install some packages and synchronize the list of
available packages with what is available on online repositories. You will then
be prompted to close the terminal. Close it and open it again. Then, again run:
```
pacman -Syu
```

This will upgrade any remaining packages.

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

The above command will install the required programs in the `/mingw64/bin`
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

Now that all the dependencies have been installed, download bmDCA.
```
git clone https://github.com/sudorook/bmDCA.git
cd bmDCA
```

Then, compile and install bmDCA in a globally accessible directory (default:
`/usr/local`) by running:
```
./autogen.sh --prefix=/usr/local && \
make -j4 && \
make install
```

Depending on your platform, the `make install` command may fail due to
permissions issues. To remedy this you can either run `sudo make install`, or
you can specify a different installation directory that does not require
administrator privileges. The latter option is particularly useful when working
on remote system not under your control.

To go with the latter option and install bmDCA in a local directory, for
example `$HOME/.local`, adjust the installation command as follows:
```
./autogen.sh --prefix=${HOME}/.local && \
make -j4 && \
make install
```

You can replace the value to the right of `--prefix=` with any other path.
Note, that you should check that it is on your system PATH.

In the event you wish to uninstall `bmDCA`, simply run `sudo make uninstall`
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
 - `-i`: training input MSA, FASTA format
 - `-I`: (_optional_) validation input MSA, FASTA format
 - `-d`: directory where output files are written
 - `-c`: (_optional_) config file for bmDCA run hyperparameters, such as
         `example/bmdca.conf`
 - `-n`: training input MSA, numerical format
 - `-N`: (_optional_) validation input MSA, numerical format
 - `-w`: training MSA sequence weights
 - `-W`: (_optional_) validation MSA sequence weights
 - `-h`: print usage information
 - `-f`: force a restart of inference loop (i.e., start at step 1)

The mapping from amino acids to integers is defined in the following way. Amino
acids are ordered as in the following string "-ACDEFGHIKLMNPQRSTVWY". They are
then mapped to the integer corresponding to their position in the string, minus
one. The gap symbol is mapped to 0, A is mapped to 1, etc...

To generate a sequence weight file, you can use the provided `process_msa` tool
bundled in bmDCA. Instructions are in a subsequent section.

__Important:__ The MSA processing function does not handle gaps represented by
'.' characters.

#### Example 1: run from FASTA file

To learn a FASTA-formatted multiple sequence alignment (with re-weighting) and
a config file:

```
bmdca -i <input_alignment.fasta> -w <input_alignment_weignts.txt> \
  -d <output_directory> -c <config_file.conf>
```

#### Example 2: restarting runs

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

**IMPORTANT** To guarantee that inferences loop produce the same results
irrespective of whether they were stopped and restarted or ran continuously,
bmDCA will check that the hyperparameters used previously and presently match.
The only fields that may be adjusted before restarting are:

1. `save_period` - the interval at which the state of the program is saved to
   disk.
2. `save_best_step` - boolean flag to save any step that produces the lowest
   RMS error with respect to the MSA so far in the inference.
3. `step_max` - the max number of steps (increase to continue a loop after it
   ends)
4. `stop_mode` - the method used to determine the convergence threshold
   (description below)
5. `stop_threshold` - the convergence threshold (lower it to continue)

The hyperparameters used by any given run are stored in the `bmdca_params.conf`
file. For boolean configuration options, it is fine to specify `1` or `0`
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

If instead, you have parameters in binary format, run:
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

### Conversion from binary to text (`arma2ascii`)

If using the `output_binary=true` flag during the inference step (recommended),
the output will be stored in an Armadillo-specific binary format. While this
allows for reproducible outputs in stopped-and-restarted inference runs, the
format is not accessible for other programs. You can use the `arma2ascii` tool
to convert binary-stored outputs to ASCII.

To convert parameters:
```
arma2ascii -p <parameters_h.bin> -P <parameters_J.bin>
```

To convert stats files:
```
arma2ascii -s <MC_stat_file.bin>
```

The output fill be stored in a `.txt` file.

### Conversion from numeric to FASTA (`numeric2fasta`)

Output sampled alignments are in numeric format. To convert them back to FASTA,
use:
```
numeric2fasta -n <numeric_msa_file> -o <output_fasta_file>
```

Presently, only numeric alignments with 21 states (20 amino acids + 1 for gaps)
can be converted to FASTA format.

### Alignment pre-preprocessing (`process_msa`)

A MSA will need to be preprocessed before any sort of inference is performed.
This is to ensure that the alignment contains mostly the relevant variation on
protein sequences. Often, for large alignment, most positions are gaps, and
most sequences consist of gaps. Such alignments will not yield informative
models.

To preprocess the MSA, you can use `process_msa` program. Command line flags
are:
 - `-i`: input multiple sequence alignment, FASTA format
 - `-n`: input multiple sequence alignment, numerical format
 - `-w`: (_optional_) sequence weight file for the input MSA
 - `-r`: (_optional_) flag to re-weight sequences before any processing is done
         (if no weights file is provided)
 - `-t`: sequence similarity threshold for re-weighting
         (if no weights file is provided)
 - `-g`: maximum allowable fraction of gaps in sequences so that highly gapped
         sequences are removed
 - `-G`: maximum allowable fraction of gaps in positions so that highly gapped
         positions are removed
 - `-s`: maximum sequence identity fraction, where sequences above the
         threshold will be pruned from the alignment
 - `-p`: positions (indexed from 0) to protect from being pruned from the
         alignment
 - `-q`: sequences (indexed form 0 in the alignment) to protect from being
         pruning
 - `-o`: output name of the processed alignment
 - `-O`: (_optional_) output name for sequence weights
 - `-h`: print usage information

Below are some examples:

To prune sequences and positions with gaps above `0.2`:
```
process_msa -i <input msa> -g .2 -G .2 -o <output file>
```

To remove sequences above 90% similarity but protect the first two sequences:
```
process_msa -i <input msa> -o <output file> -s .9 -q 0 -q 1
```

__Note:__ If pruning position and sequence gaps, gapped positions are removed
_first_ and _then_ the sequence gap proportions are computed using the pruned
set of positions.

## Configuration file options

Inference and sampling runs can be configured using a text file (see
`example/bmdca.conf`). The fields in the file are as follows:

### [bmDCA]

1. `count_max` - maximum number of iterations for Boltzmann learning process
   (default: 2000)
2. `save_period` - save parameters every `save_period` number of steps
   (default: 20)
3. `save_best_steps` - save steps that yield the lowest RMSE of the gradient
   (default: false)
4. `stop_mode` - error criterion for stopping (default: "threshold")
   * "threshold" - use the value in `error_max` below
   * "stderr" - use threshold from multinomial standard error of the MSA
   * "stderr_adj" - use threshold from multinomial standard error of the MSA
     with correction for correlated probabilities
   * "msaerr" - estimate threshold from error by bootstrapping separate subsets
     MSA against one another
5. `stop_threshold` - manual exit threshold for error convergence (default:
   1e-05)
6. `cross_validate` - subset the alignment to train one subset and validate the
   model against the other subset to assess overfitting. _If the validation
   alignment is provided at the command line, this flag is ignored and no
   subsetting is performed._
7. `validation_seqs` - number of effective sequences to sequester for
   cross-validation. _If the validation alignment is provided at the command
   line, this flag is ignored and no subsetting is performed._
8. `random_seed` - initial seed for the random number generator (default: 1)
9. `output_binary` - flag to output data in binary format, which is faster and
   more precise (default: true)
10. `update_rule` - sampler mode, 'mh' for Metropolis-Hastings and 'z-sqrt' or
    'z-barker' for Zanella, 2019. 'z-sqrt' corresponds to a balancing function
    of `sqrt(t)`, and 'z-barker' corresponds to `t/(1+t)`. (default: "mh")
11. `update_burn_times` - flag to check tune the burn times during runtime
    using MCMC sample energies and autocorrelations (walkers > 1) or lookahead
    sampling of a small set of sequences (walkers == 1). (default: true)
12. `burn_in_start` - initial burn-in time (default: 10000)
13. `burn_between_start` - initial wait time between sampling sequences
    (default: 100)
14. `adapt_up_time` - multiple to increase MCMC wait/burn-in time (default:
    1.5)
15. `adapt_down_time` - multiple to decrease MCMC wait/burn-in time (default
    0.6)
16. `step_important_max` - maximum number of importance sampling steps
    (default: 0, i.e.importance sampling disabled)
17. `coherence_min` - (default=.9999)
18. `use_ss` - flag to sample sequences equal to the effective number of
    sequences in the alignment (default: false)
19. `samples_per_walk` - number of sequences to sample for each MCMC replicate
    (default: 1000)
20. `walkers` - number of independent MCMC replicates (default: 10)

#### [[original]]

This model is a plain Boltzmann-machine without any model reparametrization and
individual learning rates for each fields and coupling. This is the default:

1. `lambda_reg_h` - regularization strength for fields (default: 0.01)
2. `lambda_reg_J` - regularization strength for couplings (default: 0.01)
3. `alpha_reg` - relative weighting of L1 (`alpha_reg=0`) and L2
   (`alpha_reg=1`) regularization (default: 1)
4. `initial_params` - choice of initializing the parameters (`profile` for an
   independent site profile model, zero otherwise)  (default: 'profile')
5. `set_zero_gauge` - set zero (Ising) gauge for parameters (default: false)
6. `allow_gap_couplings` - allow inference of gap couplings (default: true)
7. `epsilon_h` - initial values for the learning rates for fields (default:
   0.01)
8. `epsilon_J` - initial values for the learning rates for couplings (default:
   0.001)
9. `learn_rate_h_min` - smallest possible values for learning rates for fields
   (default: 0.0001) 
10. `learn_rate_h_max` - largest possible values for learning rates for fields
    (default: 0.5)
11. `learn_rate_J_min` - smallest possible values for learning rates for
    couplings (default: 0.0001) 
12. `learn_rate_J_max` - largest possible values for learning rates for
    couplings (default: 0.5)
13. `adapt_up` - scaling factor by which learning rates can be increased
    (default: 1.2)
14. `adapt_down` - scaling factor by which learning rates are decreased. Keep
    `adapt_up` x `adapt_down` less than 1 (default: 0.6)
15. `use_pos_reg` - use position-specific regularizer derived from relative
    entropy gradient for coupling frequencies in the MSA  (default: false)

#### [[reparametrization]]

Re-parametrized model described in Figliuzzi, 2018.

1. `lambda_reg_h` - regularization strength for fields (default: 0.01)
2. `lambda_reg_J` - regularization strength for couplings (default: 0.01)
3. `alpha_reg` - relative weighting of L1 (`alpha_reg=0`) and L2
   (`alpha_reg=1`) regularization (default: 1)
4. `initial_params` - choice of initializing the parameters (`profile` for an
   independent site profile model, zero otherwise)  (default: 'profile')
5. `set_zero_gauge` - set zero (Ising) gauge for parameters (default: false)
6. `epsilon_h` - initial values for the learning rates for fields (default:
   0.01)
7. `epsilon_J` - initial values for the learning rates for couplings (default:
   0.001)
8. `learn_rate_h_min` - smallest possible values for learning rates for fields
   (default: 0.0001) 
9. `learn_rate_h_max` - largest possible values for learning rates for fields
   (default: 0.5)
10. `learn_rate_J_min` - smallest possible values for learning rates for
    couplings (default: 0.0001) 
11. `learn_rate_J_max` - largest possible values for learning rates for
    couplings (default: 0.5)
12. `adapt_up` - scaling factor by which learning rates can be increased
    (default: 1.2)
13. `adapt_down` - scaling factor by which learning rates are decreased. Keep
    `adapt_up` x `adapt_down` less than 1 (default: 0.6)

#### [[adam]]

Adam (adaptive moments) is a momentum-based variant of stochastic gradient
descent where momentum terms are scaled by the second moment, so as to reduce
or enhance movement in directions that do little or enhance model accuracy.
Fast, but tends to overfit small alignments very quickly.

1. `lambda_reg_h` - regularization strength for fields (default: 0.01)
2. `lambda_reg_J` - regularization strength for couplings (default: 0.01)
3. `alpha_reg` - relative weighting of L1 (`alpha_reg=0`) and L2
   (`alpha_reg=1`) regularization (default: 1)
4. `initial_params` - choice of initializing the parameters (`profile` for an
   independent site profile model, zero otherwise)  (default: 'profile')
5. `set_zero_gauge` - set zero (Ising) gauge for parameters (default: false)
6. `allow_gap_couplings` - allow inference of gap couplings (default: true)
7. `learn_rate_h` - learning rate for fields (default: 0.01)
8. `learn_rate_J` - learning rate for couplings (default: 0.01)

#### [[adamw]]

Modified version of Adam that uses weight decay instead of L2 regularization.
The rationale is that the regularization terms in ordinary Adam get scaled up
or down by the moments and are less able to affect inference. Weight decay
(after a transform of the original L2 loss minimization equation) is
independent of the moments.

1. `lambda_decay_h` - weight decay strength for fields (default: 0.01)
2. `lambda_decay_J` - weight decay strength for couplings (default: 0.01)
3. `initial_params` - choice of initializing the parameters (`profile` for an
   independent site profile model, zero otherwise)  (default: 'profile')
4. `set_zero_gauge` - set zero (Ising) gauge for parameters (default: false)
5. `allow_gap_couplings` - allow inference of gap couplings (default: true)
6. `learn_rate_h` - learning rate for fields (default: 0.01)
7. `learn_rate_J` - learning rate for couplings (default: 0.01)
8. `eta_min` - minimum scaling value for annealing schedule (default: 0.1)
9. `eta_max` - maximum scaling value for annealing schedule (default: 1)
10. `anneal_schedule` - rule for scaling the learning rates (default: none)
    - `cos`, cosine-based annealing
    - `trap`, trapezoidal (warm-up, hot, then cool-down) annealing
11. `anneal_scale` - scaling factor for the period, extending (or shortening)
    the schedules over time and iterations (default: 2)
12. `anneal_period` - period for cycling the annealing schedule (default: 40)
13. `anneal_warm` - warm-up number of iterations to scale from `eta_min` to
    `eta_max` (default: 20)
14. `anneal_hot` - number of iterations to run at `eta_max` (default: 0)
15. `anneal_cool` - number of iterations to decrease from `eta_max` to
    `eta_min` (default: 0)

#### [[radam]]

Variant of Adam with rectified initial iterations so that the early updates,
where there are too few points to compute the second moment, follow SGDM, and
later ones use Adam's adaptive moments.

1. `lambda_reg_h` - regularization strength for fields (default: 0.01)
2. `lambda_reg_J` - regularization strength for couplings (default: 0.01)
3. `alpha_reg` - relative weighting of L1 (`alpha_reg=0`) and L2
   (`alpha_reg=1`) regularization (default: 1)
4. `initial_params` - choice of initializing the parameters (`profile` for an
   independent site profile model, zero otherwise)  (default: 'profile')
5. `set_zero_gauge` - set zero (Ising) gauge for parameters (default: false)
6. `allow_gap_couplings` - allow inference of gap couplings (default: true)
7. `learn_rate_h` - learning rate for fields (default: 0.01)
8. `learn_rate_J` - learning rate for couplings (default: 0.01)

#### [[sgdm]]

Stochastic gradient descent with momentum.

1. `lambda_reg_h` - regularization strength for fields (default: 0.01)
2. `lambda_reg_J` - regularization strength for couplings (default: 0.01)
3. `alpha_reg` - relative weighting of L1 (`alpha_reg=0`) and L2
   (`alpha_reg=1`) regularization (default: 1)
4. `initial_params` - choice of initializing the parameters (`profile` for an
   independent site profile model, zero otherwise)  (default: 'profile')
5. `set_zero_gauge` - set zero (Ising) gauge for parameters (default: false)
6. `allow_gap_couplings` - allow inference of gap couplings (default: true)
7. `beta_h` - decay rate for gradient running average for fields (default: 0.9)
8. `bata_J` - decay rate for gradient running average for couplings (default:
   0.9)
9. `learn_rate_h` - learning rate for fields (default: 0.01)
10. `learn_rate_J` - learning rate for couplings (default: 0.01)

### [bmDCA_sample]

1. `resample_max` - number of time to attempt to resample a set of decorrelated
   MCMC sequences before giving up (default: 20)
2. `random_seed` - initial seed for the random number generator (default: 1)
3. `save_interim_samples` - flag to save intermediate samples in between
   resamplings for ideal aggregate sequence properties (default: true)
4. `update_rule` - sampler mode, 'mh' for Metropolis-Hastings and 'z-sqrt' or
   'z-barker' for Zanella, 2019. 'z-sqrt' corresponds to a balancing function
   of `sqrt(t)`, and 'z-barker' corresponds to `t/(1+t)`. (default: "mh")
5. `burn_in_start` - initial burn-in time (default: 100000)
6. `burn_between_start` - initial wait time between sampling sequences
   (default: 1000)
7. `update_burn_times` - flag to check tune the burn times during runtime using
   MCMC sample energies and autocorrelations (walkers > 1) or lookahead
   sampling of a small set of sequences (walkers == 1). (default: true)
8. `adapt_up_time` - multiple to increase MCMC wait/burn-in time (default: 1.5)
9. `adapt_down_time` - multiple to decrease MCMC wait/burn-in time (default
   0.6)
10. `temperature` - temperature at which to sample sequences (default: 1.0)

## Output files

`bmdca` will output files during the course of its run:
 - `bmdca_run.log`: table of values and measurements taken at every iteration
   of the learning procedure.
 - `bmdca_params.conf`: a list of the hyperparameters used in the learning
   procedure.
 - `energies_%d.txt`: energies of each MCMC sequence, grouped by replicate
 - `samples_%d.txt`: sequences sampled from MCMC, grouped by replicate
 - `msa_numerical.txt`: numerical representation on input MSA for training
 - `msa_validate_numerical.txt`: numerical representation on input MSA for
   validation
 - `parameters_%d.txt`: learned Potts model parameters (J and h)
 - `parameters_h_%d.bin` and `parameters_J_%d.bin`: learned Potts model
   parameters (J and h), stored in arma binary format (see `output_binary=true`
   flag from config file).
 - `gradients_%d.txt`: learned Potts model gradients (J and h)
 - `gradients_h_%d.bin` and `gradients_J_%d.bin`: learned Potts model gradients
   (J and h), stored in arma binary format (see `output_binary=true` flag from
   config file).
 - `moment1_%d.txt`: first moment (J and h) of the gradients
 - `moment1_h_%d.bin` and `moment1_J_%d.bin`: first moment (J and h) of the
   gradients, stored in arma binary format (see `output_binary=true` flag from
   config file).
 - `moment2_%d.txt`: second moment (J and h) of the gradients
 - `moment2_h_%d.bin` and `moment2_J_%d.bin`: second moment (J and h) of the
   gradients, stored in arma binary format (see `output_binary=true` flag from
   config file).
 - `learning_rates_%d.txt`: per-parameter learning rates (J and h)
 - `learning_rates_h_%d.bin` and `learning_rates_J_%d.bin`: learning rates
   stored in arma binary format (see `output_binary=true` flag from config
   file).
 - `msa_weights.txt`: weights for each sequence, either a number between 0 and
   1 based on sequence similarity or 1 if re-weighting was not specified
 - `msa_validate_weights.txt`: weights for each validation sequence, either a
   number between 0 and 1 based on sequence similarity or 1 if re-weighting was
   not specified
 - `msa_stat_1p.bin`: table of frequencies for each amino acid at each position
   in the MSA
 - `msa_stat_2p.bin`: table of frequencies for pairs of amino acids at each pair
   of positions in the MSA (due to symmetry, only the 'upper triangle' of
   positions is stored)
 - `samples_stat_1p_%d.bin`: table of frequencies for each amino acid at each
   position of the set of MCMC-sampled sequences.
 - `samples_stat_1p_sigma_%d.bin`: table of standard deviation of frequencies over
   replicates for each amino acid at each position of the set of MCMC-sampled
   sequences.
 - `samples_stat_2p_%d.bin`: table of frequencies for pairs of amino acids at each
   pair of positions from the set of MCMC-sampled sequences
 - `samples_stat_2p_sigma_%d.bin`: table of standard deviation over replicates of
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

### Run log

The columns in the run log correspond to:
1. `step` - current iteration number
2. `walkers` - number of independent sampling trajectories
3. `samples-per-walk` - number of samples to take from each walker (MCMC)
4. `burn-in` - burn-in time for walkers
5. `burn-between` - burn time between samples along a trajectory
6. `total-corr` - total sequence correlation among all sampled sequences
7. `auto-corr` - correlations among sequences sampled along a trajectory
8. `cross-corr` - correlations among sequences from different trajectories
9. `auto-cross-err` - combined deviation for correlations within and between
   trajectories
10. `energy-start-avg` - mean energy after `burn-in` steps (start of
    trajectory)
11. `energy-start-sigma` - standard deviation of energies after `burn-in` steps
    (start of trajectory)
12. `energy-end-avg` - mean energy for sequences sampled at the end of
    trajectories
13. `energy-end-sigma` - standard deviation of sequences energies at the end of
    trajectories
14. `energy-err` - combined deviations for starting (after burn-in) and ending
    energies
15. `train-err-1p` - RMSE for 1p frequencies for training MSA
16. `train-err-2p` - RMSE for 2p frequencies for training MSA
17. `train-err-tot` - total RMSE for training MSA
18. `train-err-tot-min` - smallest total RMSE found so far for the training MSA
19. `validate-err-1p` - RMSE for 1p frequencies for validation MSA
20. `validate-err-2p` - RMSE for 2p frequencies for validation MSA
21. `validate-err-tot` - total RMSE for validation MSA
22. `validate-err-tot-min` - smallest total RMSE found so far for the
    validation MSA
23. `diff-avg-energy` - difference between mean training and validation MSA
    energies
24. `seed` - random seed at the current step (used for restoring states when
    re-starting bmDCA)
25. `duration` - elapsed time for the current step

Some columns may or may not be present in your log file depending on the
settings used for inference. For example, if only 1 sequence is sampled from
each trajectory (`samples_per_walk=1`), then the auto- and cross-correlation
stats along with the energy ones will be not be computed.

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
