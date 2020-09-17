

class TOOLS:

    def XY2num(self,string,X=None,Y=None):
        '''Replaces 'X' and/or 'Y' by numbers

        Replaces 'X' or 'Y' with X and Y, respectively - or both, i
        depending on which is defined.

        Parameters
        ----------
        string : str 
            string to be modified
        X : int/float/str, optional (default: None)
            value to replace 'X' with: cast to string
        Y : int/float/str, optional (default: None)
            value to replace 'Y' with: cast to string
    
        Returns
        -------
        Modified string, or input string if no replacements are
        possible
        '''


        if (X is None) and (Y is None): # No changes requested, return 
                                        # input
            return string
        elif Y is None:   # If no Y-values, only replace X
                return string.replace('$',str(X))
        elif X is None: # If no X-values, only replace y
                return string.replace('&',str(Y))
        else:   # Both X and Y are replaced
            return string.replace('$',str(X)).replace('&',str(Y))


    def file2list(self,path,fname):
        '''Parses input file to a list

        Returns a list of lines based on fname with comment and 
        empty lines removed as well as the card and subcard 
        locations in a list

        Parameters
        ----------
        path : str 
            relative directory for fname
        fname : str
            name of file to be parsed

        Returns
        -------
        ret,cards,subcards : list,list,list
            ret : str
                data input line as elements
            cards : int
                index of card headers in ret as elements
            subcards : int
                index of subcard headers in ret as elements
        '''

        ret=[]
        # Parse out comments and empty lines and store data to memory
        with open('{}/{}'.format(path,fname)) as f:  # Open the file
            for l in f:
                if len(l.split())==0: # Ignores empty lines
                    continue
                elif l.strip()[0]=='#': # Lines starting with hashes are 
                                        # considered comments and ignored
                    continue
                elif '#' in l:
                    ret.append(l.split('#')[0]) # Only store data appearing
                                                # before a hash

                else:   # Data lines are appended to list
                    ret.append(l)
    
        # Locate the cards in the list (separated by '**')
        cards=[i for i,x in enumerate(ret) if x[0:2]=='**']
        cards.append(len(ret))  # Append the end of the list for searching
                                # between entries

        # Locate the subcards in the list (separated by '*')
        subcards=[i for i,x in enumerate(ret) if x[0]=='*'] 
        subcards.append(len(ret))  # Append the end of the list for
                                   # searching between entries

        return ret,cards,subcards # Return the list of lines, the lists of
                                  # lines containing the starts of the
                                  # cards and subcards



