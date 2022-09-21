/*
 * Discount.h --
 *	Discounting schemes
 *
 * Copyright (c) 1995-2010 SRI International, 2013 Andreas Stolcke, Microsoft Corp.  All Rights Reserved.
 *
 * @(#)$Header: /home/srilm/CVS/srilm/lm/src/Discount.h,v 1.24 2019/09/09 23:13:12 stolcke Exp $
 *
 */


#ifndef _Discount_h_
#define _Discount_h_

#include "Boolean.h"
#include "File.h"
#include "Array.h"
#include "Debug.h"

#include "NgramStats.h"

/*
 * Discount --
 *	A method to manipulate counts for estimation purposes.
 *	Typically, a count > 0 is adjusted downwards to free up
 *	probability mass for unseen (count == 0) events.
 *  代码阅读参考链接：https://blog.csdn.net/kevinfight/article/details/6423930?spm=1001.2014.3001.5501
 *  Discount.h Discount.cc 这两个文件主要实现了最重要的几个折扣算法, 包括 
    a. Katz Discount (基于 Good-Turing Discounting)
    b. Absolute Discounting
    c. Natural law of succession [Eric Sven Ristad, 1995]
    d. Additive Discounting [Lidstone-Johnson-Jeffrey]
    e  Witten-Bell Discounting
    f. Kneser-Ney  Discounting
    g. Modified Kneser-Ney Discounting, [Chen, Goodman, 1998]
 *  默认情况下日志等级为0， 表示所有的信息都打印出来
 *  Discount 类是所有折扣方法的父类
 */
class Discount: public Debug {
public:
    Discount() : 
        interpolate(false) {

    };

    virtual ~Discount() {

    };

    /* discount coefficient for count */
    // 目前还不明白这个函数的意义
    virtual double discount(Count count, Count totalCount, Count observedVocab) { 
        return 1.0; 
    };	 


	/*
	 * By default, we discount float counts by discounting the
	 * integer ceiling value.
	 */
    virtual double discount(FloatCount count, FloatCount totalCount, Count observedVocab) { 
        return discount((Count)ceil(count), 
                        (Count)ceil(totalCount),
						observedVocab); 
    };

    /* weight given to the lower-order
    	* distribution when interpolating
	* high-order estimates (none by default) */
    virtual double lowerOrderWeight(Count totalCount, Count observedVocab,
					    Count min2Vocab, Count min3Vocab) { 
        return 0.0; 
    }		    

    virtual double lowerOrderWeight(FloatCount totalCount, Count observedVocab,
					                Count min2Vocab, Count min3Vocab) { 
        return lowerOrderWeight((Count)ceil(totalCount), 
                                observedVocab,
						        min2Vocab, 
                                min3Vocab); 
    };
    
    /* check if discounting disabled */
    virtual Boolean nodiscount() { 
        return false; 
    };
	
    /* save parameters to file */
    virtual void write(File &file) {

    };

    /* read parameters from file */   
    virtual Boolean read(File &file) { 
        return false; 
    };

	/*
	 * dummy estimator for when there is nothing to estimate
	 */				    
    virtual Boolean estimate(NgramStats &counts, unsigned order) { 
        return true; 
    };

	/*
	 * by default, don't allow discount estimation from fractional counts
	 */
    virtual Boolean estimate(NgramCounts<FloatCount> &counts, unsigned order) { 
        dout() << "discounting method does not support float counts\n";
        return false; 
    };

    virtual void prepareCounts(NgramCounts<NgramCount> &counts,
				unsigned order, unsigned maxOrder) { 
        return; 
    };

    virtual void prepareCounts(NgramCounts<FloatCount> &counts,
				unsigned order, unsigned maxOrder) { 
        return; 
    };

    Boolean interpolate;
    
protected:
    static unsigned vocabSize(Vocab &vocab);	/* compute effective vocabulary size */
};


/*
 * GoodTuring --
 *	The standard discounting method based on count of counts
 * 古德-图灵平滑算法是 n-gram 中标准的平滑算法，这里假设 GT 算法中最大只平滑GT_defaultMaxCount阶
 * 《自然语言处理综论》
 */

const Count GT_defaultMinCount = 1;
const Count GT_defaultMaxCount = 5;
class GoodTuring: public Discount {
public:
    GoodTuring(unsigned mincount = GT_defaultMinCount,
	           unsigned maxcount = GT_defaultMaxCount);

    double discount(Count count, Count totalCount, Count observedVocab);
    Boolean nodiscount();

    void write(File &file);
    Boolean read(File &file);

    /**
     * 平滑算法的平滑接口
     * **/
    Boolean estimate(NgramStats &counts, unsigned order);

protected:
    Count minCount;		    /* counts below this are set to 0 */        // 将低于该计数的词概率设置为0
    Count maxCount;		    /* counts above this are unchanged */       // 高于该计数的词概率不变

    Array<double> discountCoeffs;   /* cached discount coefficients */
};

/*
 * ConstDiscount --
 *	Ney's method of subtracting a constant <= 1 from all counts
 * 	(also known as "Absolute discounting").
 *	Note: this method supports interpolating higher and lower-order
 *	estimates.
 */
class ConstDiscount: public Discount
{
public:
    ConstDiscount(double d, double mincount = 0.0)
	: _discount(d < 0.0 ? 0.0 : d > 1.0 ? 1.0 : d),
	  _mincount(mincount) {};

    double discount(Count count, Count totalCount, Count observedVocab)
      { return (count <= 0) ? 1.0 : (count < _mincount) ? 0.0 : 
					(count - _discount) / count; };
    double discount(FloatCount count, FloatCount totalCount,
							Count observedVocab)
      { return (count <= 0.0) ? 1.0 :
		(count < _mincount || count < _discount) ? 0.0 : 
					(count - _discount) / count; };

    double lowerOrderWeight(Count totalCount, Count observedVocab,
					    Count min2Vocab, Count min3Vocab)
      { return _discount * observedVocab / totalCount; }

    Boolean nodiscount() { return _mincount <= 1.0 && _discount == 0.0; } ;

    Boolean estimate(NgramCounts<FloatCount> &counts, unsigned order)
      { return true; }	    /* allow fractional count discounting */

protected:
    double _discount;		    /* the discounting constant */
    double _mincount;		    /* minimum count to retain */
};

/*
 * NaturalDiscount --
 *	Ristad's natural law of succession
 */
class NaturalDiscount: public Discount
{
public:
    NaturalDiscount(double mincount = 0.0)
	: _vocabSize(0), _mincount(mincount) {};

    double discount(Count count, Count totalCount, Count observedVocab);
    Boolean nodiscount() { return false; };

    Boolean estimate(NgramStats &counts, unsigned order);
    Boolean estimate(NgramCounts<FloatCount> &counts, unsigned order) 
       { return false; };

protected:
    unsigned _vocabSize;	    /* vocabulary size */
    double _mincount;		    /* minimum count to retain */
};

/*
 * AddSmooth --
 *	Lidstone-Johnson-Jeffrey's smoothing: add a constant delta to the
 *	occurrence count of each vocabulary item.
 *
 *		p = (c + delta) / (T + N * delta) 
 *
 *	where c is the item count, T is the total count, and N is the
 *	vocabulary size. This is equivalent to a discounting factor of
 *
 *		d = (1 + delta/c) / (1 + N * delta / T)
 */
class AddSmooth: public Discount
{
public:
    AddSmooth(double delta = 1.0, double mincount = 0.0)
	: _delta(delta < 0.0 ? 0.0 : delta),
	  _mincount(mincount) {};

    double discount(Count count, Count totalCount, Count observedVocab)
      { return (count <= 0) ? 1.0 : (count < _mincount) ? 0.0 : 
					(1.0 + _delta/count) /
					(1.0 + _vocabSize*_delta/totalCount); }
    double discount(FloatCount count, FloatCount totalCount,
							Count observedVocab)
      { return (count <= 0.0) ? 1.0 : (count < _mincount) ? 0.0 : 
					(1.0 + _delta/count) /
					(1.0 + _vocabSize*_delta/totalCount); }

    Boolean nodiscount() { return _mincount <= 1.0 && _delta == 0.0; } ;

    Boolean estimate(NgramStats &counts, unsigned order)
      { _vocabSize = vocabSize(counts.vocab); return true; }
    Boolean estimate(NgramCounts<FloatCount> &counts, unsigned order)
      { _vocabSize = vocabSize(counts.vocab); return true; }

protected:
    double _delta;		    /* the additive constant */
    double _mincount;		    /* minimum count to retain */
    unsigned _vocabSize;	    /* vocabulary size */
};

/*
 * WittenBell --
 *	Witten & Bell's method of estimating the probability of an
 *	unseen event by the total number of 'new' events overserved,
 *	i.e., counting each observed word type once.
 *	Note: this method supports interpolating higher and lower-order
 *	estimates.
 */
class WittenBell: public Discount
{
public:
    WittenBell(double mincount = 0.0) : _mincount(mincount) {};

    double discount(Count count, Count totalCount, Count observedVocab)
      { return (count <= 0) ? 1.0 : (count < _mincount) ? 0.0 : 
      			((double)totalCount / (totalCount + observedVocab)); };
    double discount(FloatCount count, FloatCount totalCount,
							Count observedVocab)
      { return (count <= 0) ? 1.0 : (count < _mincount) ? 0.0 : 
      			((double)totalCount / (totalCount + observedVocab)); };
    double lowerOrderWeight(Count totalCount, Count observedVocab,
					    Count min2Vocab, Count min3Vocab)
      { return (double)observedVocab / (totalCount + observedVocab); };

    Boolean nodiscount() { return false; };

    Boolean estimate(NgramStats &counts, unsigned order)
      { return true; } ;
    Boolean estimate(NgramCounts<FloatCount> &counts, unsigned order)
      { return true; } ;
	
protected:
    double _mincount;		    /* minimum count to retain */
};

/*
 * KneserNey --
 *	Regular Kneser-Ney discounting
 */
class KneserNey: public Discount {
public:
    KneserNey(unsigned mincount = 0, 
	      Boolean countsAreModified = false,
	      Boolean prepareCountsAtEnd = false)
      : minCount(mincount),
       discount1(0.0),
       countsAreModified(countsAreModified),
       prepareCountsAtEnd(prepareCountsAtEnd) {
        // 默认构造函数，什么事情都没做
        cout << "minCount: " << minCount << endl
             << "countsAreModified: " << countsAreModified << endl
             << "prepareCountsAtEnd: " << prepareCountsAtEnd << endl;
    };

    virtual double discount(Count count, Count totalCount, Count observedVocab);
    virtual double lowerOrderWeight(Count totalCount, Count observedVocab,
					    Count min2Vocab, Count min3Vocab);
    virtual Boolean nodiscount() { return false; };

    virtual void write(File &file);
    virtual Boolean read(File &file);

    virtual Boolean estimate(NgramStats &counts, unsigned order);

    virtual void prepareCounts(NgramCounts<NgramCount> &counts, unsigned order,
							    unsigned maxOrder);

protected:
    Count minCount;		/* counts below this are set to 0 */

    double discount1;		/* discounting constant */

    Boolean countsAreModified;	/* low-order counts are already modified */
    Boolean prepareCountsAtEnd;	/* should we modify counts after computing D */
};


/*
 * ModKneserNey --
 *	Modified Kneser-Ney discounting (Chen & Goodman 1998)
 */
class ModKneserNey: public KneserNey
{
public:
    ModKneserNey(unsigned mincount = 0, 
		 Boolean countsAreModified = false,
		 Boolean prepareCountsAtEnd = false)
      : KneserNey(mincount, countsAreModified, prepareCountsAtEnd),
	discount2(0.0), discount3plus(0.0) {};

    double discount(Count count, Count totalCount, Count observedVocab);
    double lowerOrderWeight(Count totalCount, Count observedVocab,
					    Count min2Vocab, Count min3Vocab);
    Boolean nodiscount() { return false; };

    void write(File &file);
    Boolean read(File &file);

    Boolean estimate(NgramStats &counts, unsigned order);

protected:
    double discount2;		    /* additional discounting constants */
    double discount3plus;
};

#endif /* _Discount_h_ */

