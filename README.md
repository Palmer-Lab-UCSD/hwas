# hwas


## Installation and Compiling

As the package depends on the systems htslib
you need to tell R where to find the header and
library files.  To do this set the following 
environment variables

```
export HTSLIB_LIBS=-L<PATH_TO_LIB>
export HTSLIB_CFLAGS='-std=c++17 -Wall \
    -isystem<PATH_TO_HEADER_DIR>
```

Note here the use of `isystem` instead of `-I` to 
specify the path of the header files.  The reason 
for this is that we want the header files in R 
packages to be discovered before packages locally
installed on our system.

## AI Disclaimer

The AI Claude 4.7 Opus by Anthropic was used to review 
code, architectural recommendations / discussions, and
in very few cases contributed code.  Any code contributed
by Claude will be made known in the code comments or
in the git logs.

## Copyright notice

