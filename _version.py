#!/usr/bin/env python

import os
import subprocess as sp

def get_version():
    try:
        res = sp.check_output(['git', 'describe'], universal_newlines=True)
    except sp.CalledProcessError as e:
        # probably no tag to use as a reference
        # "fatal: No names found, cannot describe anything"
        raise(e, "Are there any tags for this repo?")

    res = res.strip()

    #
    # 0.1            <- no commits past the tag, so we're on the tagged commit
    # 0.1-1-g123456  <- one commit past the tag, on commit 123456. The "g" is for "git"
    #
    toks = res.split('-')
    if len(toks) == 1:
        version = toks[0]
        build = 0
    else:
        assert len(toks) == 3
        version = toks[0]
        build = toks[1]

    return version, build

if __name__ == "__main__":
    import sys
    version, build = get_version()

    usage = 'usage: {0} [version|build]'.format(sys.argv[0])
    try:
        if sys.argv[1] == 'build':
            print(build)
        elif sys.argv[1] == 'version':
            print(version)
        else:
            print(usage)
            sys.exit(1)
    except IndexError:
        print(usage)
