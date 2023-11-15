s/\bndshape(\([^)]*\))/\1.shape/g
s/\bnddata(/self.__class__(/g
s/\(isinstance(.*\)\bnddata\s*)/\1type(self))/g
