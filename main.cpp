#include <ctime>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cmath>
#include <iomanip>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/rtl/utils.hpp>
#include <boost/rtl/table_delta.hpp>
#include <boost/rtl/rename_delta.hpp>
#include <boost/rtl/range_join_delta.hpp>
#include <boost/rtl/selection_delta.hpp>
#include <boost/rtl/key_index_delta.hpp>
#include <boost/rtl/groupby_delta.hpp>
#include <boost/rtl/groupby_functors.hpp>
#include <boost/rtl/expression_registry.hpp>
#include <boost/rtl/indexed_set_impl.hpp>
#include <boost/rtl/merge_delta.hpp>

using namespace std;
using namespace boost;
using namespace rel;

///////////

static const double pi = 3.1415927;
static const double radius = 6371.0; // in km.

class EarthCoord 
{
public:
    EarthCoord(double lat, double lng) : lat_((lat/180.)*pi), lng_((lng/180.)*pi) {}
    EarthCoord() : lat_(0.), lng_(0.) {}

    double distance(const EarthCoord& other) const
    {
        //double s =  
        //    cos(lng_)*cos(lat_)*cos(other.lng_)*cos(other.lat_)+
        //    sin(lng_)*cos(lat_)*sin(other.lng_)*cos(other.lat_)+
        //    sin(lat_)*sin(other.lat_);

        //return radius*acos(s);

        double x = 69.1 * (other.lat_ - lat_); 
        double y = 69.1 * (other.lng_ - lng_) * cos(lat_/57.3);
        return sqrt(x * x + y * y);
    }

    double getLatitude() const {return lat_;}
    double getLongitude() const {return lng_;}
private:
    double lat_;
    double lng_;
};

///////////

BOOST_RTL_DEFINE_COLUMN(int, number);
BOOST_RTL_DEFINE_COLUMN(double, latitude);
BOOST_RTL_DEFINE_COLUMN(double, longitude);
BOOST_RTL_DEFINE_COLUMN(string, name);
BOOST_RTL_DEFINE_COLUMN(string, city);
//BOOST_RTL_DEFINE_COLUMN(string, state);
//BOOST_RTL_DEFINE_COLUMN(string, country);
//BOOST_RTL_DEFINE_COLUMN(string, iso_country);
//BOOST_RTL_DEFINE_COLUMN(string, continent);

struct myinfo : table_info<
    mpl::vector5<number, latitude, longitude, name, city/*, state, country, iso_country, continent*/>,
    mpl::vector1<number>
>{};

typedef table<myinfo> mytable;
//struct mytable : table<myinfo> {};

typedef mytable::value_type mytuple;

struct a;
typedef alias<number, a> number_a;
typedef alias<latitude, a> latitude_a;
typedef alias<longitude, a> longitude_a;
typedef alias<name, a> name_a;
typedef alias<city, a> city_a;

//////////

struct distance_less_than
{
    distance_less_than(double d) : d_(d) {}

    template <class It> bool operator()(const It& it) const
    {
        if (it[number()] == it[number_a()])
            return false;

        EarthCoord ec1(it[latitude()], it[longitude()]);
        EarthCoord ec2(it[latitude_a()], it[longitude_a()]);

        return ec1.distance(ec2) < d_;
    }           
    double d_;
};

//////////

class find_latitude_at_distance
{
public:
    find_latitude_at_distance(double d)
        : d_(d * pi / 20000)
    {}
	template<class It, class Table>
		typename Table::const_iterator operator()(const It& it, const Table& t) const
	{
        row<mpl::vector1<latitude_a> > sub(it[latitude()] + d_);
		return t.lower_bound(sub);
	}
private:
    double d_;
};

//////////

template<class T>
void load_table(T& t, istream& in)
{
    cout << "loading table... ";
    cout.flush();

    int cnt = 0;

    while (!in.eof())
    {
        try
        {
            string line;
            getline(in, line);

            if (line[0] == '#' || line.empty())
                continue;

            mytuple tp;

            tokenizer<char_separator<char> > tokens(line, char_separator<char>("\t"));
            tokenizer<char_separator<char> >::iterator it = tokens.begin();
            tp[number()] = lexical_cast<int>(*it++);
            tp[latitude()] = lexical_cast<double>(*it++);
            tp[longitude()] = lexical_cast<double>(*it++);
            tp[name()] = *it++;
            tp[city()] = *it++;
            //tp[state()] = *it++;
            //tp[country()] = *it++;
            //tp[iso_country()] = *it++;
            //tp[continent()] = *it;
            
            t.insert(tp);

            //if (++cnt == 50) 
            //    break;
        }
        catch (const bad_lexical_cast&)
        {}
    }
    t.begin();
    cout << count(t) << " records loaded" << endl;
}

//////////

// this counter will be sorted in the DESC order
class desc_counter_type
{
public:
    desc_counter_type(int cnt = 0) 
        : cnt_(cnt)
    {}
    void operator++()
    {
        ++cnt_;
    }
    void operator--()
    {
        --cnt_;
    }
    operator int() const
    {
        return cnt_; 
    }
    bool operator<(const desc_counter_type& other) const
    {
        return other.cnt_ < cnt_; // here
    }
private:
    int cnt_;
};

// this is a custom groupby functor
typedef column_name<counter_t<desc_counter_type> > desc_counter;

//////////

template<class T> 
void print_result(const T& t)
{
    BOOST_AUTO(begin, t.begin());
    BOOST_AUTO(end, begin);

    for (int i = 0; i < 50; ++i)
        ++end;

    for (BOOST_AUTO(it, begin); it != end; ++it)
    {
        cout << it[number()] << '\t';
        cout << setw(20) << it[city_a()] << '\t';
        cout << setw(25) << it[name_a()] << '\t';
        cout << it[desc_counter()] << '\t';
        cout << endl;
    }
}

//////////

main()
{
    ifstream in("aslatlong.txt");
    assert(in.is_open());

    mytable t;
    load_table(t, in);

    clock_t t0 = clock();

    // index the table on latitude -- bring in number for uniqueness
    BOOST_AUTO(t_by_lat, (
        materialized_index<mpl::vector2<latitude_a, number_a> >(auto_rename<a>(t)))
        );

    // self-join the table removing pairs that are _obviously_ farther away than 5 km
    BOOST_AUTO(rj,
        range_join(t, t_by_lat, find_latitude_at_distance(-5), find_latitude_at_distance(5))
        );

    // do the proper selection
    BOOST_AUTO(sel,
        selection(rj, distance_less_than(5))
        );

    // get the number of ASs
    BOOST_AUTO(gby,
        (groupby<1, mpl::vector1<desc_counter> >(sel))
        );

    // add back the AS information 
    BOOST_AUTO(mg,
        merge<1>(gby, auto_rename<a>(t))
        );

    // sort on number of ASs & city -- bring in the number to preserve uniqueness
    BOOST_AUTO(idx, (
        materialized_index<mpl::vector3<desc_counter, city_a, number> >(mg, indexed_set_strategy())
        ));

    // sort/print result
    print_result(idx);

    clock_t t1 = clock();
    cout << "The query took " << (t1-t0)*1000/CLOCKS_PER_SEC << " ms" << endl;

    // update
    transaction tr;
    expression_registry exprs;
    exprs.add(idx);

    // remove the first record in the result
    tr.remove(t, *t.lower_bound(row<mpl::vector1<number> >(3317)));
    tr.commit(exprs);

    print_result(idx);

    clock_t t2 = clock();
    cout << "Update and query re-run took " << (t2-t1)*1000/CLOCKS_PER_SEC << " ms" << endl;

    // find a record with fewer neighbours and remove it
    BOOST_AUTO(it,
        (idx.upper_bound(row<mpl::vector1<desc_counter> >(10)))
        );

    cout << it[desc_counter()] << endl;

    tr.remove(t, *t.lower_bound(row<mpl::vector1<number> >(it[number()])));
    tr.commit();

    // run the query again
    print_result(idx);

    clock_t t3 = clock();
    cout << "Update and query re-run took " << (t3-t2)*1000/CLOCKS_PER_SEC << " ms" << endl;

    return 0;
}
