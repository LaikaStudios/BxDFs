# Installation

## Requirements
* `make`
* `python3` (used by the [cpp/Makefile](../cpp/Makefile) to generate the final plugin name)
* [Pixar's RenderMan](https://renderman.pixar.com)

## Quick Start
These instructions are based on a [Rocky Linux](https://rockylinux.org/) system. Modify accordingly for your specific OS.

1. Obtain a [Professional or Non-Commercial license](https://renderman.pixar.com/store) for RenderMan to give you download access.

1. [Install RenderManProServer](https://rmanwiki.pixar.com/space/REN/542212198/Installation+and+Licensing). Default location is `/opt/pixar`.

1. [Download the RenderMan Examples](https://renderman.pixar.com/forum/downloads) and unpack these into the RenderManProServer installation location.

    For example (replace the version number and file names and locations with those matching the RenderMan version you're installing):
    ```bash
    cd /opt/pixar/RenderManProServer-27.0
    sudo tar -xzvf /home/$USER/Downloads/PixarRenderMan-Examples-27.0_2386582-linuxRHEL9_gcc11icx232.x86_64.tgz
    ```

1. Set these environment variables appropriately. These are required by the [make](https://www.gnu.org/software/make/manual/) system that's used to compile and install the plugins:
    * PIXAR_ROOT
    * RMAN_VERSION

    For example, if your version of RenderManProServer is installed in
    `/opt/pixar/RenderManProServer-27.0`, then using `bash` shell:

    ```bash
    export PIXAR_ROOT="/opt/pixar"
    export RMAN_VERSION="27.0"
    ```
    
    Since RenderManProServer requires an RMANTREE environment variable to be set to its installation location, you can conveniently use these to define it as well:
    
    ```bash
    export RMANTREE="${PIXAR_ROOT}/RenderManProServer-${RMAN_VERSION}"
    ```

1. Download or [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) this repository.
1. `cd` into the dowloaded or cloned repository's directory.

    At this point, you can use the `make` or `make all` command (they are equivalent) to build the plugins.
    You can also `cd cpp` and `make` the plugins there.

    `make clean` and `make help` can also be executed from either the top-level directory or the cpp directory.
    `make help` provides additional information about the make system and how it's controlled.

1. In order for [RenderMan](https://rmanwiki.pixar.com/display/REN) and a [RenderMan Bridge Application](https://renderman.pixar.com/bridge-tools) to find the built plugins, set this [RenderMan environment variable](https://rmanwiki.pixar.com/space/REN/542237140/Environment+Variables) appropriately:

    - RMAN_RIXPLUGINPATH

    For example, if you downloaded or cloned this repository to `${HOME}/BxDFs`, then using `bash` shell:

    ```bash
    export RMAN_RIXPLUGINPATH="${HOME}/BxDFs/build/${RMAN_VERSION}/plugins:${RMAN_RIXPLUGINPATH}"
    ```
