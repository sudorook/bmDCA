# Installation

Dependencies (installation instructions detailed below):

- [GCC](https://gcc.gnu.org/) that supports the C++11 standard and
  [OpenMP](https://en.wikipedia.org/wiki/OpenMP)
- AutoTools
- pkg-config
- [Armadillo](http://arma.sourceforge.net/)

## Dependencies

**GCC** is used to compile the source code (and dependencies, if necessary). The
code relies on the `fopenmp` flag for parallelization, so GCC is preferred over
Clang. It also needs support for the C++17 standard, so any GCC later than
version 7.4 will suffice.

**AutoTools** are a set of programs used to generate makefiles for
cross-platform compilation and installation.

**pkg-config** is a program that provides a simple interface between installed
programs (e.g. libraries and header files) and the compiler. It's used by
AutoTools to check for dependencies before compilation.

**Armadillo** is a C++ linear algebra library. It's used for storing data in
matrix structures and performing quick computations in the bmDCA inference loop.
To install, again look to your package repository.

### Linux

To install the dependencies in Linux, simply use your distributions package
manager. Commands for Debian/Ubuntu and Arch Linux are provided below:

#### Debian/Ubuntu

Run:

```sh
sudo apt-get update
sudo apt-get install git gcc g++ automake autoconf pkg-config \
  libarmadillo-dev libopenblas-dev libarpack++2-dev
```

#### Arch Linux

For Arch Linux, GCC and the build tools should have been installed with the
`base` and `base-devel` metapackages (`sudo pacman -S base base-devel`), but if
not installed, run:

```sh
sudo pacman -S gcc automake autoconf pkgconf
```

Next, install the `git` and `superlu` packages from the repos.

```sh
sudo pacman -S git superlu
```

The final dependency, Armadillo, is not in the Arch package repositories. You
can build the package from source (not recommended on rolling release distros)
or use the AUR package. Run:

```sh
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

```sh
xcode-select --install
```

This may already be installed.

Next, install Homebrew. From the [online instructions](https://brew.sh), run:

```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

If you run into permissions errors when installing Homebrew, complaining that
root owns the `/usr/local/` directory, you can change the ownership by running:

```sh
sudo chown -R <user> /usr/local/
```

where `<user>` should be substituted with your username, e.g. `john`. Your
username should be apparent form the command prompt, but if unsure, run
`whoami`.

Once Homebrew is installed, run:

```sh
brew install gcc automake autoconf pkg-config armadillo
```

This will install the most recent GCC (11.1.0 as of writing) along with
AutoTools and pkg-config.

**IMPORTANT:** The default `gcc`, located in `/usr/bin/gcc` is actually aliased
to `clang`, which is another compiler. While in principle Clang can be used to
compile bmDCA, the default version of Clang is not compatible with the `fopenmp`
compiler flag that is used to enable parallelization. Additionally, libraries
installed via Homebrew are not by default known to `pkg-config` or the linker.

Addressing issues involves overriding the `CC` and `CXX` environmental variables
to use the Homebrew-installed GCC, updating `PKG_CONFIG_PATH` with paths to any
relevant \*.pc files, and updating `LD_LIBRARY_PATH` with the paths of all the
shared object libraries linked at compile time.

Doing this for the first time is a bit bewildering, so for convenience, use the
`rcparams` file in the `tools/` directory in this repository. In it are a few
helper functions (`pkgconfig_find()` and `ld_lib_add()`) to automate the process
of finding where build dependencies are installed on your system.

Simply append the contents of that file to your shell run commands. If you don't
know what shell you're using, run:

```sh
echo $SHELL
```

For Bash, copy the contents of `rcparams` to `${HOME}/.bashrc`, and for Zsh,
copy to `${HOME}/.zshrc`. macOS versions <=10.14 (Mojave and earlier), uses Bash
as the default shell, and for >=10.15 (Catalina and later), Apple switched the
default shell to Zsh.

To copy the `rcparams` contents:

1. Append the `rcparams` file by copy-pasting the code from your favorite text
   editor, or
2. Run (as appropriate):

```sh
cat tools/rcparams >> ${HOME}/.bashrc
```

or

```sh
cat tools/rcparams >> ${HOME}/.zshrc
```

**Note:** Run commands are only executed when the shell starts, _not_ when the
files are edited. To update your currently-running shell to reflect new changes,
you can either run in the command prompt:

```sh
source ${HOME}/.bashrc
```

or simply open a new terminal. (For remote systems, you can just log out and log
in again.)

**Bash users:** Your macOS installation may not actually source the `.bashrc`
file by default. Check that the functions are actually being loaded by running:

```sh
LC_ALL=C type pkgconfig_find
```

If the function is not found, check that the `${HOME}/.bash_profile` file
exists. In it, there should be a line that looks like:

```sh
[ -f $HOME/.bashrc ] && . $HOME/.bashrc
```

If no such line is there or if no such file exists, add it and open a new
terminal.

### Windows

Before starting, install [MSYS2](https://www.msys2.org). This program is a
package distribution for GNU/Unix tools that can be used to build programs for
Windows.

The installer defaults work fine, and if prompted, open the "MSYS2" shell in the
dialog window.

Once MSYS2 is installed and open, update the base libraries by running:

```sh
pacman -Syu
```

This will download and install some packages and synchronize the list of
available packages with what is available on online repositories. You will then
be prompted to close the terminal. Close it and open it again. Then, again run:

```sh
pacman -Syu
```

This will upgrade any remaining packages.

Next, install the dependencies for bmDCA:

```sh
pacman -S git nano vim \
  autoconf automake-wrapper pkg-config make \
  mingw-w64-x86_64-toolchain \
  mingw-w64-x86_64-openmp \
  mingw-w64-x86_64-arpack \
  mingw-w64-x86_64-lapack \
  mingw-w64-x86_64-openblas \
  mingw-w64-x86_64-armadillo
```

The above command will install the required programs in the `/mingw64/bin`
directory. Unfortunately, this directory is not on the default PATH in Windows.
You will need to add it manually.

Open your `.bashrc` file in a text editor (e.g. `vim ~/.bashrc`). Nano and Vim
(two different text editors) were installed in the above command block.

Once open, add the line (at the end of the file):

```sh
export PATH="/mingw64/bin:$PATH"
```

Then, close and open the MSYS2 terminal again.

_Optionally, edit the `/etc/pacman.conf` file. Uncomment the line `#Color` and
add a new line `ILoveCandy`. Just a cosmetic flourish for `pacman`._

## Installing bmDCA (all platforms)

Now that all the dependencies have been installed on your operating systems,
download bmDCA.

```sh
git clone https://github.com/sudorook/bmDCA.git
cd bmDCA
```

Then, compile and install bmDCA in a globally accessible directory (default:
`/usr/local`) by running:

```sh
./autogen.sh --prefix=/usr/local && \
  make -j$(nproc) && \
  make install
```

Depending on your platform, the `make install` command may fail due to
permissions issues. To remedy this you can either run `sudo make install` to
install bmDCA system-wide, or you can specify a different installation directory
that does not require administrator privileges. The latter option is
particularly useful when working on remote systems not under your control.

To go with the latter option and install bmDCA in a local directory, for example
the `$HOME/.local` directory, adjust the installation command as follows:

```sh
./autogen.sh --prefix=${HOME}/.local && \
  make -j$(nproc) && \
  make install
```

You can replace the value to the right of `--prefix=` with any other directory.
Note, that you should check that it is on your system PATH.

In the event you wish to uninstall `bmDCA`, simply run `sudo make uninstall` or
`make uninstall` as appropriate.

## Run bmDCA

Test the installation by running in the terminal:

```sh
bmdca
```

If the installation worked correctly, this will print the usage information,
e.g.:

```txt
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

For full details and usage examples, see the [README](./README.md).

### Additional notes for Windows 10

Though `bmdca` should now be easily invoked from within MSYS2, one can update
the system PATH variable to make the binaries accessible system-wide, such as
from the command prompt or other terminal emulators. To update the PATH:

1. Type 'env' in the start search bar.
2. Click 'Edit the system environment variables'.
3. Click on 'Environment Variables...' toward the bottom of the window that
   opens.
4. Select 'Path' in one of the two selection windows (either 'User variables' or
   'System variables' is fine)
5. Once 'Path' is highlighted, click 'Edit...'
6. Enter the `/usr/local/bin` as a new PATH entry. You can either:
   - Click 'New' in the new window and enter the path to `/usr/local/bin` in the
     MSYS2 installation folder (default: `C:\msys64\usr\local\bin`).
   - Click the 'Browse...' button and navigate to the `C:\msys64\usr\local\bin`
     directory.
7. When the new entry is added, click 'OK' on all the opened windows to set all
   the changes. You will need to close and re-open terminals for the changes to
   be reflected.
