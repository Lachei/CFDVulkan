// Definition of sample test
#include "catch_config.hpp"
#include "catch.hpp"

int Factorial( int number ) {
    return number <= 1 ? number : Factorial( number - 2 ) * number;  // fail
    // return number <= 1 ? 1      : Factorial( number - 1 ) * number;  // pass
}

TEST_CASE( "Factorials of 1 and higher are computed (pass)", "[factorial test]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
}
