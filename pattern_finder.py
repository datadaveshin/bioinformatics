#!/Users/dave/anaconda/bin/python

import re  # regular expression module
import string
import sys
import getopt
import operator # For sorting dictionary entries by value using .itemgetter()


PATTERN_CYSTEINE = ("cysteine", re.compile("C"))
PATTERN_NGLYCOSYLATION = ("n_glycosylation", re.compile("N[^P][ST][^P]"))


class Clustal_alignment_sequence:
    """
    Makes sequence objects from supplied header and clustal alignment sequence.
    The length is also taken, but note, this is the actual length of the
    protein. Remember, the alignment sequence may contain dashes from inserts.
    """

    def __init__(self, header, aln_sequence, length=0):
        """
        Intializes the object with header, alignment sequence and length of
        sequence reported in the clustal alignment (the actual or filtered
        sequence)

        Then takes in the object's alignment sequence, and removes all the 
        dashes that are inserted into the sequence during alignment if there are 
        gaps, and sets that as the objects self.seq attibute.

        During the process two dictionaries are made to convert back and forth
        from alignment (original) sequence to a normal (filtered) sequence.
        """
        # Initialize basic header, aln_seq and length attributes
        self.header = header
        self.aln_seq = aln_sequence
        self.aln_seq_length = len(self.aln_seq)

        # Convert to a normal sequence and set as attibute 
        # Make conversion dictionaries 
        filtered_seq = ""
        orig_to_filt_dict = {}
        filt_to_orig_dict = {}
        orig_pos = 0
        filt_pos = 0
        for character in aln_sequence:
            if character in string.ascii_letters:
                orig_to_filt_dict[orig_pos] = filt_pos
                filt_to_orig_dict[filt_pos] = orig_pos
                orig_pos += 1
                filt_pos += 1
                filtered_seq += character
            elif character == "-":
                orig_pos += 1
            else: 
                print ("Error:", character, "not allowed in sequence", 
                       "\n", 
                       "Sequence should contain '-' and alphabetical",
                       "characters") 
        self.seq = filtered_seq
        self.seq_length = len(self.seq)
        if int(length) == 0:
            print ("Error: No sequence length from Clustal alignment was given")
        if self.seq_length != int(length):
            print ("Error: The calculated sequence length does not match",
                   "clustal reported length")
        self.orig_to_filt_dict = orig_to_filt_dict
        self.filt_to_orig_dict = filt_to_orig_dict
    
    def __str__(self):
        """
        Prints basic initial Alignment_sequence class information
        """
        return ("\n" + "header: " + self.header 
              + "\n" + "alignment sequence: " + self.aln_seq
              + "\n" + "length: " + str(self.length) 
              + "\n")

    def print_obj_attributes(self):
        """
        Prints additional parameters for trouble shooting
        """
        #print self.__dict__   # this prints all attributes
        for key in self.__dict__:
            #print key
            print key, "\n", self.__dict__[key], "\n"
        # print "header:\n", self.header, "\n"
        # print "length:\n", self.length, "\n"
        # print "alignment sequence:\n", self.aln_seq, "\n"
        # print "filtered sequence:\n", self.seq, "\n"
        # print "orig_to_filt_dict:\n", self.orig_to_filt_dict, "\n" 
        # print "filt_to_orig_dict:\n", self.filt_to_orig_dict, "\n"
        # pass

    def filt_to_orig(self, contracted_pos_list):
        # print "expanded_match_list", expanded_match_list
        converted_pos_list = []
        for pos_tup in contracted_pos_list:
            converted_pos_list.append((pos_tup[0], self.filt_to_orig_dict[pos_tup[1]]))
        #print "converted_position_list", converted_position_list
        return converted_pos_list


    def orig_to_filt(self, expanded_match_list):
        # print "expanded_match_list", expanded_match_list
        converted_pos_list = []
        for pos_tup in match_pos_list:
            converted_pos_list.append((pos_tup[0], self.orig_to_filt_dict[pos_tup[1]]))
        #print "converted_position_list", converted_position_list
        return converted_pos_list
    
    def get_ind_positions(self, match_tuple):
        return_list = []
        position = int(match_tuple[1])
        for letter in match_tuple[0]:
            return_list.append((letter, position))
            #seq_match_pos_list.append((letter, position))
            position += 1
        return return_list

    def build_match_string(self, aln_str, replacement_list):
        """
        Build's a string that matches the length of the original alignment
        sequence. The string begins as a series of dashes (or other characters
        if you modify the code), and then dashes are replaced with the pattern
        matches, which are defined as single amino acids and positions 
        in a list (of lists if there is more than one). Then we return the 
        string (though we could just add it directly to the object)
        """

        # Flatten the alignment position list with list comprehension:
        flat_pos_list = [tup for sub_list in replacement_list for tup in sub_list]
        
        # Replace dashes with amino acids using thier associated positions
        for pattern_tup in flat_pos_list:
            pre_pos = int(pattern_tup[1])
            post_pos = int(pattern_tup[1] + 1)
            # This string concatenation method is supposedly the quickest way
            # to "alter" a string, string types are immutable.
            aln_str = aln_str[:pre_pos] + pattern_tup[0] + aln_str[post_pos:]
        # print aln_str
        # print self.aln_seq 
        return aln_str

    def map_pattern(self, pattern):
        """
        Find all matches in the filtered seq and generates class attributes
        that includes a final match string that aligns against the alignment 
        sequence. 
        
        First the pattern is searched against the filtered sequence so that 
        matches can be found across gaps (in an alignment sequence, gaps are
        dashes).

        During the process a pattern match, which contains the pattern match
        string and the length of the pattern, is put into a list.  Each pattern 
        is then broken down to individual amino acids with positions that are 
        put into a new sequence position match list. The list is then converted
        to another list - a new alignment position match list - so that the 
        positions can be used to map the pattern matches to the alignment 
        sequence. This is what is used to make the final output string.
        """
        # Make pattern header attribute:
        self.match_header = pattern[0]

        self.match_regex = pattern[1]

        # Find the patterns within the filtered sequence:
        match_obj = re.finditer(pattern[1], self.seq)

        # End with a list of matches called match_list:
        # Make list of matches with starting positions from the match object:
        self.match_list = [(match.group(), match.start()) for match in match_obj]
        # Expanded each match so each amino acid has a position assigned to it:
        self.seq_match_pos_list = [self.get_ind_positions(match) for match in self.match_list]
        # Convert the amino acid/position list to match the aln seq:
        self.aln_match_pos_list = [self.filt_to_orig(item) for item in self.seq_match_pos_list]

        # Generate the pattern alignment string that is mapped to the alignment str
        # Build the string of dashes as input and send
        dashes = "-" * self.aln_seq_length 
        self.match_aln_str =  self.build_match_string(dashes, 
                                                      self.aln_match_pos_list)

def read_file(input_file):
    """
    Reads in a file, returns contents one line at a time.
    """
    with open(input_file) as f:
        content = f.readlines()
        return content


def gen_seq_obj_lst(file_contents):
    """ 
    Takes the individual lines from a clustal .aln file
    Remove the first few lines.
    Make a sequence object list from the remaining.
    Sequence object will contain the sequence name (header), sequence and length
    """
    sequence_list = []
    for line in file_contents:
        # Checks to see if top line containing "ClUSTAL" or a blank line
        if (line.split(' ', 1)[0] == "CLUSTAL" 
            or line[0] == ("\n") 
            or line[0] == (" ")):
                pass
        else:
            header = line.split(None, 2)[0]
            alignment_sequence = line.split(None, 2)[1]
            length = line.split(None, 2)[2]
            aln_seq_obj = Clustal_alignment_sequence(header,
                                                     alignment_sequence,
                                                     length)
            sequence_list.append(aln_seq_obj)
    return sequence_list

def screen_obj_list(parent_aln_obj, aln_obj_lst):
    # Iterate through each obj in list, including the parent object
    # Print parent_obj.match_aln_str
    # Print parent_obj.aln_seq
    for aln_obj in aln_obj_lst:
        potential_match_list = []
        for match in parent_aln_obj.aln_match_pos_list:
            potn_match_lst = []
            potn_match_str = ""
            for residue in match:
                potn_match = (aln_obj.aln_seq[residue[1]], residue[1])
                potn_match_str += aln_obj.aln_seq[residue[1]]
                potn_match_lst.append(potn_match)

            # Check if match
            if re.match(parent_aln_obj.match_regex, potn_match_str):
                potn_match_tup = (potn_match_str, potn_match_lst, True)
            else:
                potn_match_tup = (potn_match_str, potn_match_lst, False)
            potential_match_list.append(potn_match_tup)
        aln_obj.potential_match_list = potential_match_list


def output_screening_results(parent_aln_obj, aln_obj_lst, style):
    """
    Function outputs results according to user input style from the -s flag
    Currently there are 3 versions. 
    1) changes output alignment sequence to lower case, then overwrites to
       upper case if there is a match
    2) alignment sequence is all dashes unless a site does exist
    3) alignment sequence is all dashes unless a site does not exist

    Current output is to stdout
    """
    # Print top level strings:
    print parent_aln_obj.match_aln_str
    print parent_aln_obj.aln_seq

    # Lower Case Version:
    if style == "1":
        for aln_obj in aln_obj_lst:
            if aln_obj != parent_aln_obj:
                aln_seq_lower = aln_obj.aln_seq.lower()
                build_list = []
                for potn_match_tup in aln_obj.potential_match_list:
                    if potn_match_tup[2] == True:
                        build_list.append(potn_match_tup[1])
                if build_list == []:
                    output_str = aln_seq_lower
                else:
                    output_str = aln_obj.build_match_string(aln_seq_lower,
                                                            build_list)
                print output_str
 
    # Dash Version IS a site:
    elif style == "2":
        # if aln_obj != parent_aln_obj:
        for aln_obj in aln_obj_lst:
            dashes = "-" * aln_obj.aln_seq_length
            build_list = []
            for potn_match_tup in aln_obj.potential_match_list:
                if potn_match_tup[2] == True:
                    build_list.append(potn_match_tup[1])
            if build_list == []:
                output_str = dashes
            else:
                output_str = aln_obj.build_match_string(dashes, build_list)
            print output_str

    # Dash Version NOT a site:
    elif style == "3":
        # if aln_obj != parent_aln_obj:
        for aln_obj in aln_obj_lst:
            dashes = "-" * aln_obj.aln_seq_length
            build_list = []
            for potn_match_tup in aln_obj.potential_match_list:
                if potn_match_tup[2] == False:
                    build_list.append(potn_match_tup[1])
            if build_list == []:
                output_str = dashes
            else:
                output_str = aln_obj.build_match_string(dashes, build_list)
            print output_str

    # Warning that output style does not exist:
    else:
        print 'Warning output style should be "1", "2", or "3"'


def result_stats(parent_aln_obj, aln_obj_lst):
    """
    Tallies the matches and outputs statistics
    """
    match_count = 0
    for match in parent_aln_obj.aln_match_pos_list:
        # Defining variables
        # Defining the first position of each match in the parent for output.
        parent_match_first_pos = match[0][1]
        true_match_dict = {}
        sorted_true_matches = []
        false_match_dict = {}
        sorted_false_matches = []

        total_match = 0
        true_match = 0
        false_match = 0
        # Iterate thru each alignment sequence in alignment to extract stats.
        for aln_obj in aln_obj_lst:

            # Eliminate parent object from being counted in stats.
            if aln_obj != parent_aln_obj:
                match_t_or_f = aln_obj.potential_match_list[match_count][2]
            # Stat 1) Get number of true matches vs false:
                if match_t_or_f:
                    true_match += 1
                else: 
                    false_match += 1
                total_match = true_match + false_match

            # Stat 2) Make true match dictionary and false match dictionary 
            #         with number of occurances of each key:
                if match_t_or_f:
                    true_match_str = aln_obj.potential_match_list[match_count][0]
                    if true_match_str in true_match_dict:
                        true_match_dict[true_match_str] += 1
                    else:
                        true_match_dict[true_match_str] = 1
                    # Sort the dictionary results by value into a list
                    # http://stackoverflow.com/questions/613183/
                    sorted_true_matches = sorted(true_match_dict.items(), 
                                                  key=operator.itemgetter(1),
                                                  reverse = True)
                else:
                    false_match_str = aln_obj.potential_match_list[match_count][0]
                    if false_match_str in false_match_dict:
                        false_match_dict[false_match_str] += 1
                    else:
                        false_match_dict[false_match_str] = 1
                    # Sort the dictionary results by value into a list
                    # http://stackoverflow.com/questions/613183/
                    sorted_false_matches = sorted(false_match_dict.items(), 
                                                  key=operator.itemgetter(1),
                                                  reverse = True)
  


        # Stat 1) Output stats of number of true vs false matches:
        print "For match", match_count + 1, \
               parent_aln_obj.match_list[match_count][0], "at position", \
               parent_match_first_pos, "there were:"
        print  str(true_match) + "/" + str(true_match+false_match), "or", \
               str(float(true_match)/(total_match)*100) + "%", \
              "positive hits"
        print  str(false_match) + "/" + str(true_match+false_match), "or", \
               str(float(false_match)/(total_match)*100) + "%", \
              "negative hits"
        print ""

        # Stat 2) 
        print "Postive hits:"
        for item in sorted_true_matches:
            print item[0] + ":", item[1], "occurences"
        print ""
        print "Negative hits:"
        for item in sorted_false_matches:
            print item[0] + ":", item[1], "occurences" 
        print ""
        # print "false_match_dict", false_match_dict  
        # Update the match count.
        match_count += 1


def run_it(input_file, input_pattern, output_style):
    """
    docstring
    """
    # Read input file, convert each sequence in line to a sequence object,
    # return objects in list
    contents_of_file = read_file(input_file)
    seq_obj_lst = gen_seq_obj_lst(contents_of_file)

    # Send first sequence obj in list to pattern match
    top_seq_obj = seq_obj_lst[0]
    top_seq_obj.map_pattern(input_pattern)

    # top_seq_obj.print_obj_attributes()
    # seq_obj_lst[-1].print_obj_attributes()

    # Screen others for hits
    screen_obj_list(top_seq_obj, seq_obj_lst)
    # seq_obj_lst[-1].print_obj_attributes()

    # Output results
    output_screening_results(top_seq_obj, seq_obj_lst, output_style)

    # Output stats
    result_stats(top_seq_obj, seq_obj_lst)

def help_message():
    print ('\nUsage: protein_find_sites.py -i <input_file> -p <pattern> '
                   '-s <output_style> -o <output_file>\n')
    print "Pattern options are:"
    print '"cys" for cysteine'
    print '"nglyc" for n-glycosylation \n'

def main(argv):
    input_file = ''
    search_pattern = ''
    output_style = ''
    output_file = ''
    if len(argv) == 0:
        help_message()
    try:
        opts, args = getopt.getopt(argv, "hi:p:s:o:", ["ifile=", "ipattern",
                                                       "ostyle=", "ofile="])
    except getopt.GetoptError:
        help_message()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_message()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-p", "--ipattern"):
            search_pattern = arg
            if search_pattern == "c":
                input_pattern = PATTERN_CYSTEINE
            elif search_pattern == "nglyc":
                input_pattern = PATTERN_NGLYCOSYLATION
            else:
                print "Pattern options are:"
                print '"cys" for cysteine'
                print '"nglyc" for n-glycosylation'
        elif opt in ("-s", "--ostyle"):
            output_style = arg
            run_it(input_file, input_pattern, output_style)
        elif opt in ("-o", "--ofile"):
            output_file = arg
    # print 'Input file is:', input_file
    # print 'Output file is:', output_file


if __name__ == "__main__":
    main(sys.argv[1:])
