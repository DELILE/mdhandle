# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Basic logging operations.

Produces output which is suited to interactive operation at the console.

"""

# TODO: Add additional logging targets other than, or in addition to sys.stdout
# TODO: Use Python stdlib logging module
# TODO: Prepend time and context to logged messages.

import sys
import os


class Logger(object):
    """
    Logging class for :mod:`mdhandle`.

    Parameters
    -----------
    verbose : boolean
        If ``True``, output from :class:`Logger` is verbose.

    Attributes
    -----------
    verbose : boolean
        If ``True``, output from :class:`Logger` is verbose.

    Methods
    -------
    procedure_banner(msg):
        Outputs header for long running operation.
    completion_msg(msg):
        Used at end of long running operation.
    user_message(msg):
        Basic information message.
    request(msg, input_type='raw'):
        Takes user input with variable input function.
    bad_input(msg):
        Outputs information about bad user input.
    deprecation_warning(msg):
        Outputs deprecation warning message.
    warning(msg):
        Outputs warning message.  Problem is not sufficient to throw exception.
    error(msg):
        Outputs error message and raises an :exc:`Exception`.


    """

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.len_wrapper = 80

        # TODO: Add functinality to change output depending on verbose

    def procedure_banner(self, msg):
        """
        Prints banner to announce the start of program action.

        Parameters
        ----------
        msg : string

        """
        # TODO: Method to smartly breakup text from msg into 80 char wide blobs
        print('\n' + '#'*self.len_wrapper)
        print('\t %s' % msg)
        print('#'*self.len_wrapper)

    def completion_msg(self, msg, speak=False):
        """
        Standard format for completion message.

        Parameters
        ----------
        msg : string
        speak : boolean
            If ``True``, calls system to say message.
            [default=False]

        """
        if speak is True and sys.platform == 'darwin':
            # Platform dependent: OS X only
            os.system('say "%s"' % msg)
        elif sys.platform == 'darwin':
            pass
        else:
            # Other platforms
            print('Speak option is not available outside OS X')
            print('Completion message will be printed but not spoken.')

        self.procedure_banner(msg)

    def user_message(self, msg):
        """
         Standard method for simple user message.

         Parameters
         ----------
         msg : string
            Message visible to the user.
        """
        print('+ %s' % msg)

    def request(self, msg='', input_type='raw', speak=False):
        """
        Used for printing user requests and questions.
        
        Utilizes either ``input()`` or ``raw_input()``
        functions from standard library depending on value of ``input_type``

        Parameters
        ----------
        msg : string, [default='']
            String providing context for request to user.
        input_type : string, {'raw', 'input'}, [default='raw']
           Selects the type of function used for catching user input.
        speak : boolean
            If ``True``, calls system to speak request.
            [default=False]

        """
        if speak is True and sys.platform == 'darwin':
            # Platform dependent: OS X only
            os.system('say "%s"' % msg)
        elif sys.platform == 'darwin':
            pass
        else:
            # Other platforms
            print('Speak option is not available outside OS X')
            print('Completion message will be printed but not spoken.')
        
        if input_type == 'raw':
            return raw_input('>> %s ' % msg)
        elif input_type == 'input':
            return input('>> %s ' % msg)
        else:
            self.error('@ Input type %s does not exist' % input_type)

    def bad_input(self, msg):
        """
        Printed in response to bad user input.

        Parameters
        ----------
        msg : string

        """
        print('!! Bad Input: %s' % msg)

    def warning(self, msg):
        """
        Used for printing user warning messages.

        Typically used in case of important user feed back, but when not
        so severe that exception must be thrown.

        Parameters
        ----------
        msg : string

        """
        print('!! Warning:')
        print('    %s' % msg)

    def deprication_warning(self, msg):
        """
        Used for printing deprication warnings.
        
        Does not throw any exceptions.

        Parameters
        ----------
        msg : string

        """
        # Add deprication decorator like:
        # http://code.activestate.com/recipes/391367-deprecated/
        print('!! Deprecation Warning:')
        print('    %s' % msg)

    def error(self, msg):
        """
        Print error message.

        Parameters
        ----------
        msg: string

        """
        print('!! ERROR:')
        print('    %s' % msg)


#=============================================================================

def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
