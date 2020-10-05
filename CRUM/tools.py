

class Tools:

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



    def pickplot(self,ax,pl):
        if pl=='plot':
            return ax.plot
        elif pl=='semilogy':
            return ax.semilogy
        elif pl=='semilogx':
            return ax.semilogx
        elif pl=='loglog':
            return ax.loglog
        else:
            print('Unknown plot option "{}"'.format(pl))

    def getdecade(self,x):
        from numpy import log10
        exp=max(int(log10(max(x.max(),1e-20))),
                int(log10(max(1e-20,abs(x).max()))))
        
            
        return x/(10**exp),exp

    def setkwargs(self,ax,x,y,tot,**kwargs):
        
        if 'xlabel' in kwargs: ax.set_xlabel(kwargs.get('xlabel'))
        else: ax.set_xlabel('Time [ms]')

        if 'ylim' in kwargs: ax.set_ylim(kwargs.get('ylim'))  
        else: ax.set_ylim(1.1*min(0,y.min()),1.1*max(max(tot.max(),y.max()),
                            0.1*abs(y.min())))

        if 'xlim' in kwargs: ax.set_xlim(kwargs.get('xlim'))
        else: ax.set_xlim(x[0],x[-1])

        

    def plotax( self,ax,x,y,plot='plot',nuclei=True,linewidth=2,plottot=True,
                ylabel=r'[$ 10^{{{}}}$]',totlabel='Total',
                plotmax=1e3,**kwargs):
        ''' Plots '''
        from numpy import sum

        y,exp=self.getdecade(y)
        tot=nuclei*self.totparticles(y)+sum(y,axis=0)*(not nuclei)

        if plottot is True:
            self.pickplot(ax,plot)(x,tot,'k',linewidth=linewidth,label=totlabel)

        for i in range(int(min(plotmax,len(y)))):
            self.pickplot(ax,plot)(x,y[i,:],linewidth=linewidth,label=self.slist[i])


        ax.set_ylabel(ylabel.format(exp))
        self.setkwargs(ax,x,y,tot*plottot,**kwargs)



    def plotparterror(  self,ax,x,y,plot='plot',linewidth=2,color='r',
                        ylabel=r'[$ 10^{{{}}}$]',right=False,**kwargs):
        
        y=self.totparticles(y)
        y,exp=self.getdecade(y-y[0])
        if right is True:
            ax=ax.twinx()
        self.pickplot(ax,plot)(x,y,color=color)

        ax.set_ylabel(ylabel.format(exp))
        self.setkwargs(ax,x,y,y*0,**kwargs)


    def plotEerror( self,ax,x,y,plot='plot',linewidth=2,color='r',tscal=1e3,
                    ylabel=r'[$ 10^{{{}}}$]',**kwargs):
        from numpy import sum
        
        y=sum(y,axis=0)
        y,exp=self.getdecade(y)
        self.pickplot(ax,plot)(x,y,color=color)

        ax.set_ylabel(ylabel.format(exp))
        self.setkwargs(ax,x,y,y*0,**kwargs)

    def plottotpart(self,ax,x):
        ax.axhline(self.totparticles(self.getdecade(x)[0]),color='k',linewidth=0.5)

    def setgrid(self,fig):
        for ax in fig.get_axes():
            ax.grid('g:',linewidth=0.2)


