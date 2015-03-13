# SybilLite install guide #

## Quick setup ##

You know what you're doing.  Or hate to read.  Have at it.

```
    $ mkdir /var/www/sybillite
    $ ./configure --prefix /var/www/sybillite
    $ ./install
```

Then edit your apache conf as directed by the configure script and restart apache.  Direct your browser to:

> http://localhost/sybillite


## Detailed setup ##

SybilLite is a web interface to view comparative genomics data.  To start, you'll need a web server, presumably Apache.  The CGI scripts have a set of required modules:

  * Config::IniFiles
  * Data::Dumper
  * DBD::SQLite
  * File::Basename
  * File::Copy
  * GD
  * GD::SVG
  * HTML::Template
  * JSON
  * URI::Escape

These are checked by the configuration script, so don't worry about checking them individually on your own.

For this guide we'll assume you're installing in a subdirectory of /var/www/

Start by creating the directory you want to use:

```
    $ mkdir /var/www/sybillite
```

Then run the configure script, passing that directory using the --prefix option.

```
    $ ./configure --prefix /var/www/sybillite
```

It will check that you have all the required perl modules and, if so, create an 'install' executable.  It will also print out an example entry you'll need to add to your apache.conf file.  This is because SybilLite mixes HTML and CGI documents and you need to direct Apache to always execute the CGI ones.  For our example install, that looks like this:

```
    <Directory "/var/www/sybillite">
        Options +ExecCGI
        AllowOverride Limit
        DirectoryIndex index.html index.cgi
        AddHandler cgi-script .cgi

        <FilesMatch ".ini">
                Deny from all
        </FilesMatch>
    </Directory>
```

Note that the precise file you need to put that in depends heavily on your Unix/Linux distribution in use.  For Ubuntu 11.10, for example, it goes here:

```
    /etc/apache2/sites-available/default
```

Check the documentation for your OS and apache distribution.  Once  you've added the entry, you'll need to restart Apache.

Now, you're ready to run the 'install' script created during your configuration.  It will place the files necessary to get the interface up and running.

Point your browser to:

```
    http://localhost/sybillite/
```

If you need help:

```
    http://groups.google.com/group/sybillite-users?pli=1
```

If you encounter any problems, please file an issue here:

```
    http://code.google.com/p/sybillite/issues/list
```






