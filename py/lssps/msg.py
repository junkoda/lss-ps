import numbers
import lssps._lssps as c


class LogLevel:
    def __init__(self):
        self._str = ('debug', 'verbose', 'info', 'warn', 'error', 'fatal',
                     'silent')
        self._dict = {}
        for i, s in enumerate(self._str):
            self._dict[s] = i

        self.loglevel = 'warn'

    @property
    def loglevel(self):
        n = c._msg_get_loglevel()
        return self._str(n)

    @loglevel.setter
    def loglevel(self, value):
        if isinstance(value, str):
            if not value in self._dict:
                raise TypeError('loglevel must be one of the following: ' + ' '.join(self._str))
                
            nlevel = self._dict[value]
        elif isinstance(loglevel, numbers.Integral):
            nlevel = value
        else:
            raise TypeError('loglevel must be a string or integer')

        c._msg_set_loglevel(nlevel)
    
    

