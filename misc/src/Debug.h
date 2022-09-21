/*
 * Debug.h --
 *	General object debugging facility
 *
 * Debug is a Mix-in class that provides some simple, but consistent
 * debugging output handling.
 *
 * Copyright (c) 1995-2010 SRI International.  All Rights Reserved.
 *
 * @(#)$Header: /home/srilm/CVS/srilm/misc/src/Debug.h,v 1.8 2013/03/30 15:55:25 stolcke Exp $
 *
 */

#ifndef _Debug_h_
#define _Debug_h_

#ifdef PRE_ISO_CXX
# include <iostream.h>
#else
# include <iostream>
using namespace std;
#endif

#include <Boolean.h>

/*
 * Here is the typical usage for this mixin class.
 * First, include it in the parents of some class FOO
 *
 * class FOO: public OTHER_PARENT, public FOO { ... }
 *
 * Inside FOO's methods use code such as
 *
 *	if (debug(3)) {
 *	   dout() << "I'm feeling sick today\n";
 *	}
 *
 * Finally, use that code, after setting the debugging level
 * of the object and/or redirecting the debugging output.
 *
 *      FOO foo;
 *	foo.debugme(4); foo.dout(cout);
 *
 * Debugging can also be set globally (to affect all objects of
 * all classes.
 *
 *	foo.debugall(1);
 *
 */

/**
 * SRILM 中日志类
 * 
 * ***/
class Debug {
public:
    Debug(unsigned level = 0) : 
        nodebug(false), 
        debugLevel(level), 
        debugStream(&cerr) {

    };
    virtual ~Debug() {};	    /* prevent warning about no virtual dtor */

    // 传入日志等级
    Boolean debug(unsigned level) const {/* true if debugging */
	    return (!nodebug && (debugAll >= level || debugLevel >= level)); 
    };

    /* set object's debugging level */
    // 设置当前实例的日志等级
    virtual void debugme(unsigned level) { 
        debugLevel = level; 
    };
				    
    /* set global debugging level */
    // 设置全局日志等级
    void debugall(unsigned level) { 
        debugAll = level; 
    };
				    
    unsigned debuglevel() const { 
        return debugLevel; 
    };

    /* output stream for use with << */
    virtual ostream &dout() const { 
        return *debugStream; 
    };
				    
    /* redirect debugging output */
    virtual ostream &dout( ostream &stream ) { 
        debugStream = &stream; return stream; 
    };

    Boolean nodebug;		    /* temporarily disable debugging */
private:
    static unsigned debugAll;	    /* global debugging level */
    unsigned debugLevel;	    /* level of output -- the higher the more*/
    ostream *debugStream;	    /* current debug output stream */
};

#endif /* _Debug_h_ */

