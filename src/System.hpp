#pragma once

#include <vector>
#include <boost/operators.hpp>
#include <ostream>

template< class T , size_t Dim >
class point :
    boost::additive1< point< T , Dim > ,
    boost::additive2< point< T , Dim  > , T ,
    boost::multiplicative2< point< T , Dim > , T
    > > >
    {
    public:

        const static size_t dim = Dim;
        typedef T value_type;
        typedef point< value_type , dim > point_type;

        // ...
        // constructors
        //<-
        point( void )
        {
            for( size_t i=0 ; i<dim ; ++i ) m_val[i] = 0.0;
        }

        point( value_type val )
        {
            for( size_t i=0 ; i<dim ; ++i ) m_val[i] = val;
        }

        point( value_type x , value_type y , value_type z = 0.0 )
        {
            if( dim > 0 ) m_val[0] = x;
            if( dim > 1 ) m_val[1] = y;
            if( dim > 2 ) m_val[2] = z;
        }
        //->

        // ...
        // operators
        //<-
        T operator[]( size_t i ) const { return m_val[i]; }
        T& operator[]( size_t i ) { return m_val[i]; }

        point_type& operator+=( const point_type& p )
        {
            for( size_t i=0 ; i<dim ; ++i )
                m_val[i] += p[i];
            return *this;
        }

        point_type& operator-=( const point_type& p )
        {
            for( size_t i=0 ; i<dim ; ++i )
                m_val[i] -= p[i];
            return *this;
        }

        point_type& operator+=( const value_type& val )
        {
            for( size_t i=0 ; i<dim ; ++i )
                m_val[i] += val;
            return *this;
        }

        point_type& operator-=( const value_type& val )
        {
            for( size_t i=0 ; i<dim ; ++i )
                m_val[i] -= val;
            return *this;
        }

        point_type& operator*=( const value_type &val )
        {
            for( size_t i=0 ; i<dim ; ++i )
                m_val[i] *= val;
            return *this;
        }

        point_type& operator/=( const value_type &val )
        {
            for( size_t i=0 ; i<dim ; ++i )
                m_val[i] /= val;
            return *this;
        }

        //->

    private:

        T m_val[dim];
    };

    //...
    // more operators
    //]

    //
    // the - operator
    //
    template< class T , size_t Dim >
    point< T , Dim > operator-( const point< T , Dim > &p )
    {
        point< T , Dim > tmp;
        for( size_t i=0 ; i<Dim ; ++i ) tmp[i] = -p[i];
        return tmp;
    }

    //
    // scalar product
    //
    template< class T , size_t Dim >
    T scalar_prod( const point< T , Dim > &p1 , const point< T , Dim > &p2 )
    {
        T tmp = 0.0;
        for( size_t i=0 ; i<Dim ; ++i ) tmp += p1[i] * p2[i];
        return tmp;
    }

    //
    // norm
    //
    template< class T , size_t Dim >
    T norm( const point< T , Dim > &p1 )
    {
        return scalar_prod( p1 , p1 );
    }

    //
    // absolute value
    //
    template< class T , size_t Dim >
    T abs( const point< T , Dim > &p1 )
    {
        return sqrt( norm( p1 ) );
    }

    //
    // output operator
    //
    template< class T , size_t Dim >
    std::ostream& operator<<( std::ostream &out , const point< T , Dim > &p )
    {
        if( Dim > 0 ) out << p[0];
        for( size_t i=1 ; i<Dim ; ++i ) out << " " << p[i];
        return out;
    }

typedef point< double , 3 > point_type;
typedef std::vector< point_type > container_type;
typedef std::vector< double > scalar_type;
typedef std::vector< bool > mask_type;

namespace doorstep
{
   point_type center_of_mass( const container_type &x , const scalar_type &m )
   {
       const size_t n = x.size();
       double overall_mass = 0.0;
       point_type mean( 0.0 );
       for( size_t i=0 ; i<n ; ++i )
       {
	   overall_mass += m[i];
	   mean += m[i] * x[i];
       }
       if( !x.empty() ) mean /= overall_mass;
       return mean;
   }
}
