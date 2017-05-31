# stogie
Validate Deletion SVs with CIGAR strings

## Install

#### Clone from GitHub

```
git clone --recursive https://github.com/dantaki/stogie.git
```

#### Compile with [CMake](https://cmake.org/)

```
cd stogie/
mkdir build && cd build/

cmake .. && make 
```

##### Binary executable found under `stogie/build/src/stogie`

## Usage 

Stogie currently supports validation of deletion SV

`stogie --help`

```
8""""8    ""8""    8"""88    8""""8    8     8""""
8           8      8    8    8    "    8     8
8eeeee      8e     8    8    8e        8e    8eeee
    88      88     8    8    88  ee    88    88
e   88      88     8    8    88   8    88    88
8eee88      88     8eeee8    88eee8    88    88eee

stogie         Validate Deletions with CIGAR strings
Version: 1.0	Author: Danny Antaki <dantaki@ucsd.edu>

Usage: stogie -i <in.bam> -r <sv.bed> -x <FLOAT> -q <INT> -o <output.txt>

Options:
    -i        Input: BAM filename
    -r        Input: SV bed file
    -x        Minimum reciprocal overlap [0.5]
    -q        Mapping quality threshold [10]
    -o        Output: filename
```
## Input

* `-i    Input: BAM filename`
   *  BAM file must be indexed with `samtools index`

* `-r    Input: SV bed file`
   * The first four tab-delimited columns must be `CHROM START END DEL`   

## Output

| CHR | START | END | LENGTH | OVERLAP | READNAME | ALIGNMENT | STRANDS | SV | TYPE |
| --- | ----- | --- | ------ | ------- | -------- | --------- | ------- | --- | --- | 
| 11 | 26679952 | 26680054 | 103 | 0.932039 | STOGIE_READ | + | 11:26679948-26680047 | DEL |

## Acknowledgements

stogie uses [BamTools](https://github.com/pezmaster31/bamtools)

> *BamTools*
> Copyright (c) 2009-2010 Derek Barnett, Erik Garrison, Gabor Marth, Michael Stromberg

Notice: BamTools requires [zlib](http://zlib.net/).

## License
stogie: validate Deletions with CIGAR strings 

Copyright (C) 2017 Danny Antaki

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## Contact

dantaki@ucsd.edu
