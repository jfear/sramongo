"""Some small helpers for dealing with xml."""


def valid_path(func):
    """Validates that a XML path exists.

    A decorator function that makes sure a XML path (i.e.
    xml.etree.ElementTree.ElemenTree.Element) exists. If the path does not
    exists then return an empty dictionary.

    """
    def new_func(*args, **kwargs):
        # If the current path is present
        if args[1] is None:
            print('Not valid path.', func, args)
            return {}
        else:
            return func(*args, **kwargs)

    return new_func
