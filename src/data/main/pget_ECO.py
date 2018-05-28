#!/usr/bin/env python
#
# Copyright (c) 2010 Giorgos Keramidas.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

import os
import sys

children = {}
maxjobs = 50                # maximum number of concurrent jobs
jobs = []                   # current list of queued jobs

# Default wget options to use for downloading each URL
wget = ["wget", "-q", "-nd", "-np", "-c", "-r"]

# Spawn a new child from jobs[] and record it in children{} using
# its PID as a key.
def spawn(cmd, *args):
    argv = [cmd] + list(args)
    pid = None
    try:
        pid = os.spawnlp(os.P_NOWAIT, cmd, *argv)
        children[pid] = {'pid': pid, 'cmd': argv}
    except Exception as inst:
        print("'%s': %s" % ("\x20".join(argv), str(inst)))
    print( "spawned pid %d of nproc=%d njobs=%d for '%s'" % \
        (pid, len(children), len(jobs), "\x20".join(argv)))
    return pid

if __name__ == "__main__":
    # Build a list of wget jobs, one for each URL in our input file(s).
    for fname in sys.argv[1:]:
        try:
            with open(fname,'r') as file_f:
                file_f_lines = file_f.readlines()
                for u in file_f_lines:
                    cmd = wget + [u.strip('\r\n')]
                    jobs.append(cmd)
        except IOError:
            pass
    print( "%d wget jobs queued" % len(jobs))

    # Spawn at most maxjobs in parallel.
    while len(jobs) > 0 and len(children) < maxjobs:
        cmd = jobs[0]
        if spawn(*cmd):
            del jobs[0]
    print( "%d jobs spawned" % len(children))

    # Watch for dying children and keep spawning new jobs while
    # we have them, in an effort to keep <= maxjobs children active.
    while len(jobs) > 0 or len(children):
        (pid, status) = os.wait()
        print( "pid %d exited. status=%d, nproc=%d, njobs=%d, cmd=%s" % \
            (pid, status, len(children) - 1, len(jobs), \
             "\x20".join(children[pid]['cmd'])))
        del children[pid]
        if len(children) < maxjobs and len(jobs):
            cmd = jobs[0]
            if spawn(*cmd):
                del jobs[0]
