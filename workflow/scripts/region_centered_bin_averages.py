#!/usr/bin/env python
#----------------------------------------------------------------------------------------------------------------------
#
# Create bins around input features
#
# Created by:
#     Wanding Zhou (with documentation by Jacob Morrison)
#
# Notes:
#     Jan 2022
#         - Initial creation from https://github.com/zwdzwd/wzmetagene
#
#----------------------------------------------------------------------------------------------------------------------
import argparse
import numpy as np

# A terminology note:
#
# The "start" of the chromosome is base 0 (or 1 in 1-based indexing) in the reference. The "end" of the chromosome is,
# therefore, the last (highest numbered) base in the reference.

class Record(object):

    def __init__(self, fields, args):
        """Initialize record.

        Inputs -
            fields - list of fields from BED file
            args   - argparse.parse_args() from running script
        """
        self.fields = fields       # list of fields from BED file entry
        self.chrm = fields[0]      # chromosome
        self.beg  = int(fields[1]) # start of region
        self.end  = int(fields[2]) # end of region

        # Collapse initial region to the middle of the region
        if args.collapse:
            self.mid = (self.beg + self.end + 1) / 2 # midpoint of initial region
            self.beg = self.mid # set start to midpoint
            self.end = self.mid # set end to midpoint

        self.step1 = args.flankstep # step size in reverse direction (relative to reference and strand)
        self.step2 = args.flankstep # step size in forward direction (relative to reference and strand)

    def __str__(self):
        """String method for print()."""
        return f'chrm: {self.chrm}, start: {self.beg}, end: {self.end}, step1: {self.step1}, step2: {self.step2}'

    def sample_forward(self, args, index_func):
        """Move towards the end of the chromosome (relative to the reference) from the target region.

        Inputs -
            args       - argparse.parse_args() from running script
            index_func - function for how to adjust the index value for bin
        Returns -
            None
        """
        # If forward step is negative, then don't create any windows towards the end of the chromosome
        if self.step2 < 0:
            return

        # When moving towards the end of the chromosome, the windows begin at the 3' end (relative to the reference) of
        # the region of interest.
        # NB: this variable is a temporary variable and may be adjusted depending on provided inputs
        _window_beg = self.end

        # Create args.flanknumber of windows while moving towards the end of the chromosome
        for i in range(args.flanknumber):
            # This is a temporary variable and may be adjusted depending on inputs
            _window_end = _window_beg + self.step2

            if args.outer:
                # If args.outer, set the end of the window to be the 3' end of the temporary window:
                # (_window_end-1, _window_end)
                window_end = int(_window_end)
                window_beg = window_end - 1
            elif args.middle:
                # If args.middle, set the end of the window to be the middle of the temporary window:
                # (middle of window - 1, middle of window)
                window_mid = int((_window_beg + _window_end)/2.0)
                window_beg = window_mid-1
                window_end = window_mid
            else:
                # Otherwise use the whole temporary window
                window_beg = int(_window_beg)
                window_end = int(_window_end)

            # Index of the probe window moving towards the end of the chromosome (relative to the reference)
            index = index_func(i)

            # The index is negative for reverse strand (-) regions
            if index < 0:
                if args.flankbygene:
                    # For args.flankbygene, set the region column as the number index of the current window
                    reg = '({:d})-({:d})'.format(-i-1, -i)
                elif args.flanktoneighbor:
                    # For args.flanktoneighbor, set the region column as the corresponding percentage covered by the
                    # current window between the current interval and the next
                    if args.fold:
                        # If args.fold, halve region percentage to account for setting the forward and backward regions
                        # indices the same
                        reg = '({:d})-({:d})%'.format(
                            int(float(-i-1) / args.flanknumber / 2 * 100),
                            int(float(-i) / args.flanknumber / 2 * 100)
                        )
                    else:
                        # Otherwise, you can leave the percentage as is
                        reg = '({:d})-({:d})%'.format(
                            int(float(-i-1) / args.flanknumber * 100),
                            int(float(-i) / args.flanknumber * 100)
                        )
                else:
                    # Otherwise, set the region column as bases covered relative the end of the interval
                    reg = '({:d})-({:d})'.format(int(self.end - window_end), int(self.end - window_beg))

            # Handle forward strand intervals
            else:
                # The descriptions are the same as the negative index case, so won't rehash them here
                if args.flankbygene:
                    reg = '{:d}-{:d}'.format(i, i+1)
                elif args.flanktoneighbor:
                    # FIXME: This should behave in a similar manner to the args.flanktoneighbor above (where there is an
                    #        if-else block for args.fold. This is how it's done for sample_backward(), but for some
                    #        reason it's not done here.
                    reg = '{:d}-{:d}%'.format(
                        int(float(i) / args.flanknumber * 100),
                        int(float(i+1) / args.flanknumber * 100)
                    )
                else:
                    reg = '{:d}-{:d}'.format(int(window_beg - self.end), int(window_end - self.end))

            # When args.collapse is set, you can have index = 0
            # To be honest, I don't know what the area value is supposed to represent...
            if index >= 0:
                area = 1
            else:
                area = -1

            # Expand probe window if desired
            if args.expansion > 0:
                # Start properly accounts for running off the start of the chromosome, but end does not account for
                # running off the end of the chromosome. Maybe this can be accounted for somehow?
                window_beg = max(window_beg - args.expansion,0)
                window_end = window_end + args.expansion

            # Potential FIXME: I think we would still want to print out this information, even if window_beg == 0.
            #                  Therefore, I think that window_beg > 0 should be changed to window_beg >= 0
            if window_beg > 0 and window_end > window_beg:
                print(
                    '{}\t{:d}\t{:d}\t{:d}\t{}\t{:d}\t{}'.format(
                        self.chrm,  # chromosome of probe window
                        window_beg, # start of probe window (0-based)
                        window_end, # end of probe window (1-based, non-inclusive)
                        index,      # 0-indexed index of the probe window
                        reg,        # region covered by probe window
                        area,       # not exactly sure what this is....
                        '\t'.join(self.fields) # original values from input region of interest
                    )
                )

            # Move current temporary end to start of next probe window
            _window_beg = _window_end

    def sample_backward(self, args, index_func):
        """Move towards the start of the chromosome (relative to the reference) from target region.

        Inputs -
            args       - argparse.parse_args() from running script
            index_func - function for how to adjust the index value for bin
        Returns -
            None
        """
        # If backward step is negative, then don't create any windows towards the start of the chromosome
        if self.step1 < 0:
            return

        # When moving towards the start of the chromosome, the windows begin at the 5' end (relative to the reference)
        # of the region of interest.
        # NB: this variable is a temporary variable and may be adjusted depending on provided inputs
        _window_end = self.beg

        # Create args.flanknumber of windows while moving towards the start of the chromosome
        for i in range(args.flanknumber):
            # This is a temporary variable and may be adjusted depending on inputs
            _window_beg = _window_end - self.step1

            if args.outer:
                # If args.outer, set the end of the window to be the 5' end of the temporary window:
                # (_window_beg, _window_beg+1)
                window_beg = int(_window_beg)
                window_end = window_beg + 1
            elif args.middle:
                # If args.middle, set the end of the window to be the middle of the temporary window:
                # (middle of window - 1, middle of window)
                window_mid = int((_window_beg + _window_end)/2.0)
                window_beg = window_mid-1
                window_end = window_mid
            else:
                # Otherwise use the whole temporary window
                window_beg = int(_window_beg)
                window_end = int(_window_end)

            # Index of the probe window moving towards the start of the chromosome (relative to the reference)
            index = index_func(i)

            # The index is negative for forward strand (+) regions
            if index < 0:
                if args.flankbygene:
                    # For args.flankbygene, set the region column as the number index of the current window
                    reg = '({:d})-({:d})'.format(-i-1, -i)
                elif args.flanktoneighbor:
                    # For args.flanktoneighbor, set the region column as the corresponding percentage covered by the
                    # current window between the current interval and the previous
                    if args.fold:
                        # If args.fold, halve region percentage to account for setting the forward and backward regions
                        # indices the same
                        reg = '({:d})-({:d})%'.format(
                            int(float(-i-1) / args.flanknumber / 2 * 100),
                            int(float(-i) / args.flanknumber / 2 * 100)
                        )
                    else:
                        # Otherwise, you can leave the percentage as is
                        reg = '({:d})-({:d})%'.format(
                            int(float(-i-1) / args.flanknumber * 100),
                            int(float(-i) / args.flanknumber * 100)
                        )
                else:
                    # Otherwise, set the region column as bases covered relative the end of the interval
                    reg = '({:d})-({:d})'.format(int(window_beg - self.beg), int(window_end - self.beg))

            # Handle reverse strand intervals
            else:
                # The descriptions are the same as the negative index case, so won't rehash them here
                if args.flankbygene:
                    reg = '{:d}-{:d}'.format(i, i+1)
                elif args.flanktoneighbor:
                    if args.fold:
                        reg = '{:d}-{:d}%'.format(
                            int(float(i) / args.flanknumber / 2 * 100),
                            int(float(i+1) / args.flanknumber / 2 * 100)
                        )
                    else:
                        reg = '{:d}-{:d}%'.format(
                            int(float(i) / args.flanknumber * 100),
                            int(float(i+1) / args.flanknumber * 100)
                        )
                else:
                    reg = '{:d}-{:d}'.format(int(self.beg - window_end), int(self.beg - window_beg))

            # When args.collapse is set, you can have index = 0
            # To be honest, I don't know what the area value is supposed to represent...
            if index >= 0:
                area = 1
            else:
                area = -1

            # Expand probe window if desired
            if args.expansion > 0:
                # Start properly accounts for running off the start of the chromosome, but end does not account for
                # running off the end of the chromosome. Maybe this can be accounted for somehow?
                window_beg = max(window_beg - args.expansion,0)
                window_end = window_end + args.expansion

            # Potential FIXME: I think we would still want to print out this information, even if window_beg == 0.
            #                  Therefore, I think that window_beg > 0 should be changed to window_beg >= 0
            if window_beg > 0 and window_end > window_beg:
                print(
                    '{}\t{:d}\t{:d}\t{:d}\t{}\t{:d}\t{}'.format(
                        self.chrm,  # chromosome of probe window
                        window_beg, # start of probe window (0-based)
                        window_end, # end of probe window (1-based, non-inclusive)
                        index,      # 0-indexed index of the probe window
                        reg,        # region covered by probe window
                        area,       # not exactly sure what this is....
                        '\t'.join(self.fields) # original values from input region of interest
                    )
                )

            # Move current temporary start to end of next probe window
            _window_end = _window_beg

    def sample_internal(self, args, index_func):
        """Move within the target region.

        Inputs -
            args       - argparse.parse_args() from running script
            index_func - function for how to adjust the index value for bin
        Returns -
            None
        """
        if args.outer:
            # sentinels is a list of locations to sample from inside the region of interest
            # NB: If args.numinternal = 1, then the location for the internal region will always be:
            #     (region start, region start + 1)
            # NB: Here, start is set to be 1-based
            sentinels = list(np.linspace(self.beg+1, self.end, args.numinternal))

            for i in range(len(sentinels)):
                window_end = int(sentinels[i]) # end of window to probe
                window_beg = window_end - 1    # start of window to probe (reset to 0-based)

                # Expand probe window if desired
                if args.expansion > 0:
                    # Start properly accounts for running off the start of the chromosome, but end does not account for
                    # running off the end of the chromosome. Maybe this can be accounted for somehow?
                    window_beg = max(window_beg - args.expansion,0)
                    window_end = window_end + args.expansion

                # Index of the internal probe window
                index = index_func(i)

                # Potential FIXME: I think we would still want to print out this information, even if window_beg == 0.
                #                  Therefore, I think that window_beg > 0 should be changed to window_beg >= 0
                if window_beg > 0 and window_end > window_beg:
                    print(
                        '{}\t{:d}\t{:d}\t{:d}\t{:d}-{:d}%\t0\t{}'.format(
                            self.chrm,  # chromosome of probe window
                            window_beg, # start of probe window (0-based)
                            window_end, # end of probe window (1-based, non-inclusive)
                            index,      # 0-indexed index of the probe window
                            int(float(index)/args.numinternal*100),   # low end percent of range of the region covered by the probe window
                            int(float(index+1)/args.numinternal*100), # high end percent of range of the region covered by the probe window
                            '\t'.join(self.fields) # original values from input region of interest
                        )
                    )

            return

        # If args.outer is not set, then break up the internal region of interest into numinternal regions
        # For example, for region (0, 100), then if numinternal = 1 the probe window will cover (0, 100). If
        # numinternal = 2, the probe windows will be (0, 50) and (50, 100)
        # NB: Everything stays in BED 0-based indexing for this part of the function
        sentinels = list(np.linspace(self.beg, self.end, args.numinternal+1))

        for i in range(len(sentinels)-1):
            window_beg = int(sentinels[i])   # start of probe window
            window_end = int(sentinels[i+1]) # end of probe window

            # If args.middle is set, we want to set the probe window as the middle position of the window. So, for the
            # previously listed examples, the (0, 100) interval would be set to (49, 50) and the (0, 50) and (50, 100)
            # intervals would become (24, 25) and (74, 75)
            if args.middle:
                window_mid = int((sentinels[i] + sentinels[i+1])/2) # find middle of the current window
                window_beg = window_mid-1 # reset the start to the middle (setting value to 0-based index)
                window_end = window_mid   # reset the end to the middle location

            # Expand probe window if desired
            # See above comments for handling the end position relative to the end of the chromosome
            if args.expansion > 0:
                window_beg = max(window_beg - args.expansion,0)
                window_end = window_end + args.expansion

            # Index of the internal probe window (0 for internal window with only 1 internal probe window)
            index = index_func(i)

            # See the above potential FIXME, as this applies here
            if window_beg > 0 and window_end > window_beg:
                print(
                    '{}\t{:d}\t{:d}\t{:d}\t{:d}-{:d}%\t0\t{}'.format(
                        self.chrm,  # chromosome of probe window
                        window_beg, # start of probe window (0-based)
                        window_end, # end of probe window (1-based, non-inclusive)
                        index,      # 0-indexed index of the probe window
                        int(float(index)/args.numinternal*100),   # low end percent of range of the region covered by the probe window
                        int(float(index+1)/args.numinternal*100), # high end percent of range of the region covered by the probe window
                        '\t'.join(self.fields) # original values from input region of interest
                    )
                )

        return

def process_record(r, r0, r2):
    """Parse target region, creating windows based on the flanking inputs given on the command line.

    Inputs -
        r  - record to process
        r0 - previous record in file
        r2 - next record in file
    Returns -
        None
    """
    if args.flankbygene:
        # When flanking by gene (or more accurately, flanking by region of interest), set the step size to be equal to
        # the size of the region of interest
        # NB: specifically, this is related to the size of the windows in the region of interest, so the number of
        #     internal windows (args.numinternal) plays into this size
        r.step1 = float(r.end - r.beg + 1) / (args.numinternal)
        r.step2 = r.step1

    elif args.flanktoneighbor:
        # When flanking by neighbor, the entire space between the end of the previous region and the start of the
        # current region will be broken up into args.flanknumber number of windows (and likewise for the space between
        # the end of the current region and the start of the next region)

        # If no previous region exists or the end of the previous region overlaps the start of the current region, then
        # don't create any windows upstream (towards the start of the chromosome relative to the reference) of the
        # current region
        if r0 is None or r.beg < r0.end:
            r.step1 = -1
        else:
            # Evenly divide up region between end of previous region and start of current region into args.flanknumber
            # number of windows
            r.step1 = float(r.beg - r0.end + 1) / (args.flanknumber)

            # Based on the description of args.fold, I'm not entirely sure why the the step size is halved when
            # args.fold is set
            if args.fold:
                r.step1 /= 2

        # If no next region exists or the start of the next region overlaps the end of the current region, then don't
        # create any windows downstream (towards the end of the chromosome relative to the reference) of the current
        # region
        if r2 is None or r.end > r2.beg:
            r.step2 = -1
        else:
            # Evenly divide up region between end of current region and start of next region into args.flanknumber
            # number of windows
            r.step2 = float(r2.beg - r.end + 1) / (args.flanknumber)

            # Based on the description of args.fold, I'm not entirely sure why the the step size is halved when
            # args.fold is set
            if args.fold:
                r.step2 /= 2

    # If strand information is not included, assume all regions are on the forward (+) strand
    # Otherwise pick off the strand information from the input file
    if args.strand is None:
        strand = '+'
    else:
        strand = r.fields[args.strand-1]

    # If args.ignoreend, then the end of the region is ignored. For forward strand regions, this is the 3' end relative
    # to the reference, so the region will be set to (start, start+1). For reverse strand regions, this is the 5' end
    # relative to the reference, so the region will be set to (end-1, end)
    if args.ignoreend:
        args.numinternal = 1
        if strand == '+':
            r.end = r.beg + 1
        else:
            r.beg = r.end - 1

    # Potential FIXME: If args.fold is included, then every line will be printed twice - once for this if-block and once
    #                  for the strand if-else block. Further, the indices between the two prints will differ, as
    #                  args.fold makes the index equivalent between the upstream and downstream windows. This could be
    #                  solved by placing a return as the last expression in the block
    if args.fold:
        r.sample_backward(args, lambda i: -i-1)
        r.sample_internal(args, lambda i: min(i, args.numinternal-1-i))
        r.sample_forward(args, lambda i: -i-1)

    if strand == '+':
        # Create probe regions for forward strand (+) intervals
        r.sample_backward(args, lambda i: -i-1)
        r.sample_internal(args, lambda i: i)
        r.sample_forward(args, lambda i: args.numinternal+i)
    else:
        # For reverse strand (-) intervals, moving upstream from the CTCF cite is actually moving downstream relative to
        # the reference (i.e., the position value gets larger)
        r.sample_forward(args, lambda i: -i-1)
        r.sample_internal(args, lambda i: args.numinternal-1-i)
        r.sample_backward(args, lambda i: args.numinternal+i)

    return

def main(args):
    """Main function to run.

    Inputs -
        args   - argparse.parse_args() from running script
    Returns -
        None
    """
    # If collapsing interval to the middle, then (re)set the number of points sampled in the internal interval
    if args.collapse:
        if args.outer: # if args.outer, then keep the collapsed middle for internal probe
            args.numinternal = 1
        else: # ignore internal interval if --outer is not set
            args.numinternal = 0

    # Loop through intervals in BED file
    r0 = None # previous record (interval), only used with args.flanktoneighbor
    r  = None # current record (interval)
    r2 = None # next record (interval), only used with args.flanktoneighbor
    for line in args.table:
        r2 = Record(line.strip('\n').split('\t'), args)

        if r is not None: # skip the first line
            if r.chrm == r2.chrm: # current and next intervals fall on the same chromosome
                process_record(r, r0, r2)
            else: # next interval falls on a different chromosome
                process_record(r, r0, None)

        r0 = r
        r  = r2

    # Handle the last record
    r2 = None
    process_record(r, r0, r2)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate meta gene')
    parser.add_argument('table', help="Input bed file", type = argparse.FileType('r'), default='-')

    # How to handle intervals
    parser.add_argument(
        '--middle',
        action = 'store_true',
        help = 'use middle point of each probe window'
    )
    parser.add_argument(
        '--outer',
        action = 'store_true',
        help = 'use outer point of each probe window (w.r.t. to the target region)'
    )
    parser.add_argument(
        '--collapse',
        action = 'store_true',
        help = 'collapse internal interval to middle - all probe windows will now be relative to this location'
    )

    # Control how the probe window sizes are determined
    # Set step size to region of interest size - highest precedence when setting probe step size
    parser.add_argument(
        '--flankbygene',
        action = 'store_true',
        help = 'allow the size of steps to vary according to the region of interest length'
    )
    # Set step size to be flanknumber equal steps to next region - second highest precedence when setting step size
    parser.add_argument(
        '--flanktoneighbor',
        action = 'store_true',
        help = 'size of step is dependent on the closest region of interest to current region'
    )
    # Create equal size steps with length args.flankstep - used when flankbygene and flanktoneighbor are not set
    parser.add_argument(
        '-f', '--flankstep',
        type = int,
        default=100,
        help = 'set step size to X bases (default 100)'
    )
    parser.add_argument(
        '-m', '--flanknumber',
        type = int,
        default=30 ,
        help = 'number of points to sample outside of the region of interest (default 30)')

    # Expand probe windows - expands window in both directions, occurs for internal and external windows
    parser.add_argument(
        '--expansion',
        type = int,
        default = 0,
        help = 'number of bases to expand window, expands window in both directions (default 0)'
    )

    # Control internal sampling
    parser.add_argument(
        '-n', '--numinternal',
        type = int,
        default = 30,
        help = 'number of points to sample in the region of interest, --middle ignores this (default 30)'
    )

    # Other inputs
    parser.add_argument(
        '--fold',
        action = 'store_true',
        help = 'use the same index for intervals on both sides of the target region (usually used when strand is irrelevant)'
    )
    parser.add_argument(
        '-s', '--strand',
        type = int,
        default = None,
        help = 'the column which contains strand information in input file - if None then ignore strand (default None)'
    )
    parser.add_argument(
        '--ignoreend',
        action ='store_true',
        help = 'ignore the end of the input interval'
    )

    parser.set_defaults(func=main)
    args = parser.parse_args()

    try:
        args.func(args)
    except IOError:
        exit
