#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from __future__ import print_function

import sys
import time
import os
import imp
import shlex
import struct
import platform
import subprocess

import tqdm

__all__ = ['Progress', 'OldProgress']

class Progress(object):
    
    """
    Before I had my custom progressbar here.
    Now it is a wrapper around the great progressbar `tqdm`.
    Old implementation moved to `OldProgress` class.
    """
    
    def __init__(self, total = None, name = "Progress",
             interval = None, percent = True, status = 'initializing',
             done = 0, init = True, unit = 'it'):
        
        self.name = name
        self.interval = max(int(total / 100), 1) if interval is None else interval
        self.total = total
        self.done = done
        self.status = status
        self.unit = unit
        self.start_time = time.time()
        self.min_update_interval = 0.1
        self.last_printed_value = 0
        
        if init:
            self.init_tqdm()
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def init_tqdm(self):
        self.tqdm = tqdm.tqdm(total = self.total,
                              desc = '%s: %s' % (self.name, self.status),
                              unit_scale = True,
                              unit = self.unit)
        self.last_updated = time.time()
    
    def step(self, step = 1, msg = None, status = 'busy', force = False):
        """
        Updates the progressbar by the desired number of steps.
        
        :param int step: Number of steps or items.
        """
        self.done += step
        
        if force or (self.done % self.interval < 1.0 and \
            time.time() - self.last_updated > self.min_update_interval):
            
            self.set_status(status)
            
            this_update = max(0, self.done - self.last_printed_value)
            
            if this_update == 0:
                self.tqdm.refresh()
                self.tqdm.fp.flush()
            else:
                self.tqdm.update(int(this_update))
                
            self.last_printed_value = self.done
            self.last_updated = time.time()
    
    def terminate(self, status = 'finished'):
        """
        Terminates the progressbar and destroys the tqdm object.
        """
        self.step(self.total - self.done, force = True, status = status)
        self.tqdm.close()
    
    def set_total(self, total):
        """
        Changes the total value of the progress bar.
        """
        self.total = total
        self.tqdm.total = total
        self.step(0)
    
    def set_done(self, done):
        """
        Sets the position of the progress bar.
        """
        self.done = done
        self.tqdm.n = self.done
        self.tqdm.last_print_n = self.done
        self.step(0)
    
    def set_status(self, status):
        """
        Changes the prefix of the progressbar.
        """
        if status != self.status:
            self.status = status
            self.tqdm.set_description(self.get_desc())
            self.tqdm.refresh()
    
    def get_desc(self):
        """
        Returns a formatted string of the description, consisted of
        the name and the status. The name supposed something constant
        within the life of the progressbar, while the status is there
        to give information about the current stage of the task.
        """
        return '%s%s%s%s' % (' ' * 8,
                             self.name,
                             ' -- ' if len(self.name) else '',
                             self.status)


class OldProgress(object):
    def __init__(self,
                 total=None,
                 name="Progress",
                 interval=3000,
                 percent=True,
                 status='initializing'):
        self.status = status
        self.name = name
        self.interval = interval
        self.set_total(total)
        self.done = 0
        self.percent = percent
        self.in_ipython()
        self.last_updated = time.time()
        self.start_time = time.time()
        self.min_update_interval = 0.0 \
            if not self.ipython or self.ipython == 'terminal' \
            else 1.0
        self.width = get_terminal_size()[0] - 1
        self.clear_line()
        if self.percent:
            _ = sys.stdout.write("\r\t:: %s: %s 0.00%%" %
                                 (self.name, '%s,' % self.status))
        else:
            self.total = int(self.total)
            _ = sys.stdout.write("\r\t:: %s: %s 0/%u" %
                                 (self.name, '%s,' % self.status, self.total))
        sys.stdout.flush()

    def step(self, step=1, msg=None, status='working on it', force=False):
        self.status = status
        self.done += step
        if force or (
                self.done % self.interval < 1.0 and
                time.time() - self.last_updated > self.min_update_interval):
            self.clear_line()
            if self.percent:
                _ = sys.stdout.write("\r\t:: %s: %s %.2f%% %s" % (
                    self.name, ''
                    if self.status is None else '%s,' % self.status,
                    float(self.done) / float(self.total) * 100.0, ''
                    if msg is None else '[%s]' % msg))
            else:
                _ = sys.stdout.write("\r\t:: %s: %s %u/%u %s" % (
                    self.name, '%s,' % self.status, self.done, self.total, ''
                    if msg is None else '[%s]' % msg))
            sys.stdout.flush()
            self.last_updated = time.time()

    def set_total(self, total):
        self.total = total or 9999999999.0

    def set_done(self, done):
        self.done = done

    def set_status(self, status):
        # if self.ipython == 'notebook':
        #    status = '%s; progress indicator disabled in IPython notebook.' % status
        self.step(step=0, status=status, force=True)

    def clear_line(self):
        _ = sys.stdout.write("\r" + " " * self.width)

    def terminate(self, status='finished'):
        self.finished_time = time.time()
        self.seconds_elapsed = self.finished_time - self.start_time
        m, s = divmod(self.seconds_elapsed, 60)
        h, m = divmod(m, 60)
        self.time_elapsed = "%d:%02d:%02d" % (h, m, s)
        self.clear_line()
        self.status = status
        if self.percent:
            _ = sys.stdout.write(
                "\r\t:: %s: %s 100.0%% [%s elapsed]" %
                (self.name, '%s,' % self.status, self.time_elapsed))
        else:
            _ = sys.stdout.write("\r\t:: %s: %s %u/%u [%s elapsed]" %
                                 (self.name, '%s,' % self.status, self.total,
                                  self.total, self.time_elapsed))
        _ = sys.stdout.write("\n")
        sys.stdout.flush()
        self.last_updated = time.time()

    def in_ipython(self):
        self.ipython = False
        if 'ipykernel' in sys.modules:
            self.ipython = 'notebook'
        elif 'Ipython' in sys.modules:
            self.ipython = 'terminal'

# terminal size
# from https://gist.github.com/jtriley/1108174


def get_terminal_size():
    """ getTerminalSize()
     - get width and height of console
     - works on linux,os x,windows,cygwin(windows)
     originally retrieved from:
     http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python

    credits to: https://gist.github.com/jtriley/1108174
    """
    current_os = platform.system()
    tuple_xy = None
    if current_os == 'Windows':
        tuple_xy = _get_terminal_size_windows()
        if tuple_xy is None:
            tuple_xy = _get_terminal_size_tput()
            # needed for window's python in cygwin's xterm!
    if current_os in ['Linux', 'Darwin'] or current_os.startswith('CYGWIN'):
        tuple_xy = _get_terminal_size_linux()
    if tuple_xy is None:
        tuple_xy = (80, 25)  # default value
    return tuple_xy


def _get_terminal_size_windows():
    try:
        from ctypes import windll, create_string_buffer
        # stdin handle is -10
        # stdout handle is -11
        # stderr handle is -12
        h = windll.kernel32.GetStdHandle(-12)
        csbi = create_string_buffer(22)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
        if res:
            (bufx, bufy, curx, cury, wattr, left, top, right, bottom, maxx,
             maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
            sizex = right - left + 1
            sizey = bottom - top + 1
            return sizex, sizey
    except:
        pass


def _get_terminal_size_tput():
    # get terminal width
    # src:
    # http://stackoverflow.com/questions/263890/how-do-i-find-the-width-height-of-a-terminal-window
    try:
        cols = int(subprocess.check_call(shlex.split('tput cols')))
        rows = int(subprocess.check_call(shlex.split('tput lines')))
        return (cols, rows)
    except:
        pass


def _get_terminal_size_linux():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl
            import termios
            cr = struct.unpack('hh',
                               fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
            return cr
        except:
            pass

    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            return None
    return int(cr[1]), int(cr[0])
