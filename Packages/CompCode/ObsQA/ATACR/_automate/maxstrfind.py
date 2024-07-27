def maxstrfind(s,p):
        '''
        Wild this doesnt exist in Python libraries. Finds LAST occurence of expression in a string.
        '''
        i = 0
        while i>-1:
                io = i
                i = s.find(p,i+1,len(s))
        return io