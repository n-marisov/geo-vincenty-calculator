<?php

namespace Maris\Geo\Vincenty;

use Maris\Interfaces\Geo\Aggregate\LocationAggregateInterface;
use Maris\Interfaces\Geo\Calculator\BearingCalculatorInterface;
use Maris\Interfaces\Geo\Calculator\DistanceCalculatorInterface;
use Maris\Interfaces\Geo\Factory\LocationFactoryInterface;
use Maris\Interfaces\Geo\Finder\DestinationFinderInterface;
use Maris\Interfaces\Geo\Model\EllipsoidInterface;
use Maris\Interfaces\Geo\Model\LocationInterface;

/***
 * Калькулятор Vincenty.
 * @author Марисов Николай Андреевич
 */
class VincentyCalculator implements DistanceCalculatorInterface, BearingCalculatorInterface,DestinationFinderInterface
{
    /**
     * Пи * 2 для вычислений.
     */
    protected const M_2_PI = M_PI * 2;

    /**
     * Пи * 3 для вычислений.
     */
    protected const M_3_PI = M_PI * 3;

    /**
     * Фабрика координат.
     * @var LocationFactoryInterface
     */
    private LocationFactoryInterface $locationFactory;

    /***
     * Эллипсоид для расчета.
     * @var EllipsoidInterface
     */
    protected EllipsoidInterface $ellipsoid;

    /**
     * Максимальное число итераций.
     * @var int|null
     */
    protected ?int $iMax;

    /**
     * @param LocationFactoryInterface $locationFactory
     * @param EllipsoidInterface $ellipsoid
     * @param int|null $iMax
     */
    public function __construct( LocationFactoryInterface $locationFactory, EllipsoidInterface $ellipsoid, int $iMax = null )
    {
        $this->locationFactory = $locationFactory;
        $this->ellipsoid = $ellipsoid;
        $this->iMax = $iMax;
    }


    public function calculateDistance( LocationAggregateInterface|LocationInterface $start, LocationAggregateInterface|LocationInterface $end ): float
    {
        $start = self::deg2radLocationToArray( $start );
        $end = self::deg2radLocationToArray( $end );
        return $this->inverse( $start["lat"], $start["long"], $end["lat"], $end["long"] )["distance"];
    }

    public function findDestination(LocationAggregateInterface|LocationInterface $location, float $initialBearing, float $distance): LocationInterface
    {
        $location = self::deg2radLocationToArray( $location );
        return $this->locationFactory->new(
            ...$this->direct($location["lat"],$location["long"], deg2rad( $initialBearing ), $distance)["destination"]
        );
    }

    public function calculateInitialBearing(LocationAggregateInterface|LocationInterface $start, LocationAggregateInterface|LocationInterface $end): float
    {
        $start = self::deg2radLocationToArray( $start );
        $end = self::deg2radLocationToArray( $end );
        return $this->inverse( $start["lat"], $start["long"], $end["lat"], $end["long"] )["bearing"]["initial"];
    }

    public function calculateFinalBearing(LocationAggregateInterface|LocationInterface $start, LocationAggregateInterface|LocationInterface $end): float
    {
        $start = self::deg2radLocationToArray( $start );
        $end = self::deg2radLocationToArray( $end );
        return $this->inverse( $start["lat"], $start["long"], $end["lat"], $end["long"] )["bearing"]["final"];
    }

    /**
     * Вычисляет ряд А
     * @param float $k
     * @return float
     */
    protected function calcA( float $k ):float
    {
        return (1 +  $k ** 2 / 4) / (1-$k);
    }

    /***
     * Вычисляет ряд В
     * @param float $k
     * @return float
     */
    protected function calcB( float $k ):float
    {
        return $k * ( 1 - 3 * $k ** 2 / 8);
    }

    /**
     * Вычисляет коэффициент для расчета рядов А и В
     * @param float $uSq
     * @return float
     */
    protected function calcK( float $uSq ):float
    {
        return ( ($s = sqrt(1 + $uSq )) - 1 ) / ( $s + 1 );
    }

    /**
     * Вычисляет параметр С.
     * @param float $cosSquAlpha
     * @return float
     */
    protected function calcC( float $cosSquAlpha ):float
    {
        return $this->ellipsoid->getFlattening() / 16 * $cosSquAlpha * ( 4 + $this->ellipsoid->getFlattening()  * (4 - 3 * $cosSquAlpha) );
    }
    /**
     * Вычисляет U в квадрате.
     * @param $cosSquareAlpha
     * @return float
     */
    protected function calcUSquare( $cosSquareAlpha ):float
    {
        $squareB = $this->ellipsoid->getEquatorRadius() ** 2 ;
        return $cosSquareAlpha * ($this->ellipsoid->getPolarRadius() ** 2 - $squareB) / $squareB;
    }

    /**
     * @param float $B
     * @param float $sinSigma
     * @param float $cosSigma
     * @param float $cos2SigmaM
     * @return float
     */
    protected function calcDeltaSigma(float $B, float $sinSigma, float $cosSigma, float $cos2SigmaM):float
    {
        return $B * $sinSigma * ($cos2SigmaM + $B / 4
                * ($cosSigma * (-1 + 2 * $cos2SigmaM * $cos2SigmaM) - $B / 6
                    * $cos2SigmaM * (-3 + 4 * $sinSigma * $sinSigma)
                    * (-3 + 4 * $cos2SigmaM * $cos2SigmaM)
                )
            );
    }

    /**
     * Обратная задача
     * @param $startLat
     * @param $startLon
     * @param $endLat
     * @param $endLon
     * @return array{distance:float,bearing:array{initial:float,final:float}}
     */
    protected function inverse( $startLat, $startLon, $endLat, $endLon ):array
    {
        $L = $endLon - $startLon;

        $tanU1 = ( 1 - $this->ellipsoid->getFlattening() ) * tan($startLat);
        $cosU1 = 1 / sqrt(1 + $tanU1 * $tanU1);
        $sinU1 = $tanU1 * $cosU1;
        $tanU2 = (1 - $this->ellipsoid->getFlattening() ) * tan($endLat);
        $cosU2 = 1 / sqrt(1 + $tanU2 ** 2 );
        $sinU2 = $tanU2 * $cosU2;

        $lambda = $L;

        $iterations = 0;

        do {
            $sinLambda = sin($lambda);
            $cosLambda = cos($lambda);
            $sinSqSigma = ($cosU2 * $sinLambda) * ($cosU2 * $sinLambda)
                + ($cosU1 * $sinU2 - $sinU1 * $cosU2 * $cosLambda) * ($cosU1 * $sinU2 - $sinU1 * $cosU2 * $cosLambda);
            $sinSigma = sqrt($sinSqSigma);

            if ($sinSigma == 0) return [
                "distance" => 0,
                "bearing" => [
                    "initial" => 0,
                    "final" => 0
                ]
            ];

            $cosSigma = $sinU1 * $sinU2 + $cosU1 * $cosU2 * $cosLambda;
            $sigma = atan2($sinSigma, $cosSigma);
            $sinAlpha = $cosU1 * $cosU2 * $sinLambda / $sinSigma;
            $cosSquAlpha = 1 - $sinAlpha * $sinAlpha;

            /**
             * Устанавливаем на 0 на случай экваториальных линий
             */
            $cos2SigmaM = ($cosSquAlpha !== 0.0) ? $cosSigma - 2 * $sinU1 * $sinU2 / $cosSquAlpha : 0;


            $C = $this->calcC( $cosSquAlpha );

            $lambdaP = $lambda;
            $lambda = $L + (1 - $C) * $this->ellipsoid->getFlattening() * $sinAlpha
                * ($sigma + $C * $sinSigma * ($cos2SigmaM + $C * $cosSigma * (-1 + 2 * $cos2SigmaM * $cos2SigmaM)));
            $iterations++;

            /**
             * Выходим из цикла если превышено максимальное число итераций.
             */
            if(isset($this->iMax) && $this->iMax <=  $iterations)
                break;

        } while ( abs($lambda - $lambdaP) > 1E-12 );

        $uSq = $this->calcUSquare( $cosSquAlpha );
        $K = $this->calcK( $uSq );
        $A = $this->calcA( $K );
        $B = $this->calcB( $K );

        $a1 = atan2($cosU2 * $sinLambda, $cosU1 * $sinU2 - $sinU1 * $cosU2 * $cosLambda);
        $a2 = atan2($cosU1 * $sinLambda, -$sinU1 * $cosU2 + $cosU1 * $sinU2 * $cosLambda);

        $a1 = fmod($a1 + self::M_2_PI, self::M_2_PI);
        $a2 = fmod($a2 + self::M_2_PI, self::M_2_PI);

        return [
            "distance" => $this->ellipsoid->getEquatorRadius() * $A * ( $sigma - $this->calcDeltaSigma( $B, $sinSigma, $cosSigma, $cos2SigmaM )),
            "bearing" =>[
                "initial" => rad2deg($a1),
                "final" => rad2deg($a2)
            ]
        ];
    }

    /**
     * Реализация прямой задачи.
     * @param float $startLat
     * @param float $startLong
     * @param float $bearing
     * @param float $distance
     * @return array
     */
    protected function direct( float $startLat, float $startLong , float $bearing, float $distance ):array
    {
        $sinAlpha1 = sin( $bearing );
        $cosAlpha1 = cos( $bearing );

        $tanU1 = ( 1 - $this->ellipsoid->getFlattening() ) * tan($startLat);
        $cosU1 = 1 / sqrt(1 + $tanU1 * $tanU1);
        $sinU1 = $tanU1 * $cosU1;
        $sigma1 = atan2($tanU1, $cosAlpha1);
        $sinAlpha = $cosU1 * $sinAlpha1;
        $cosSquAlpha = 1 - $sinAlpha * $sinAlpha;


        $K = $this->calcK( $this->calcUSquare( $cosSquAlpha ) );
        $A = $this->calcA( $K );
        $B = $this->calcB( $K );


        $sigmaS = $distance / ( $this->ellipsoid->getEquatorRadius() * $A );
        $iterations = 0;
        do{
            $cos2SigmaM = cos(2 * $sigma1 + ($sigma ?? $sigma = $sigmaS) );
            $sinSigma = sin($sigma);
            $cosSigma = cos($sigma);
            $deltaSigma = $this->calcDeltaSigma( $B,  $sinSigma,  $cosSigma,  $cos2SigmaM );
            $sigmaS = $sigma;
            $sigma = $distance / ( $this->ellipsoid->getEquatorRadius() * $A) + $deltaSigma;
            $iterations++;

            /**
             * Выходим из цикла если превышено максимальное число итераций.
             */
            if(isset($this->iMax) && $this->iMax <=  $iterations)
                break;


        }while( abs($sigma - $sigmaS) > 1E-12 );

        $tmp = $sinU1 * $sinSigma - $cosU1 * $cosSigma * $cosAlpha1;
        $lambda = atan2($sinSigma * $sinAlpha1, $cosU1 * $cosSigma - $sinU1 * $sinSigma * $cosAlpha1);
        $C = $this->calcC( $cosSquAlpha );
        $L = $lambda
            - (1 - $C) * $this->ellipsoid->getFlattening() * $sinAlpha
            * ($sigma + $C * $sinSigma * ($cos2SigmaM + $C * $cosSigma * (-1 + 2 * $cos2SigmaM ** 2)));

        return [
            "destination" =>[
                rad2deg( atan2(
                    $sinU1 * $cosSigma + $cosU1 * $sinSigma * $cosAlpha1,
                    (1 - $this->ellipsoid->getFlattening()) * sqrt($sinAlpha * $sinAlpha + $tmp * $tmp)
                ) ),
                rad2deg( fmod($startLong + $L + self::M_3_PI, self::M_2_PI ) - M_PI )
            ],
            "finalBearing" => rad2deg( fmod(atan2( $sinAlpha, -$tmp ) + self::M_2_PI , self::M_2_PI ) )
        ];

    }

    /**
     * Приводит точка-подобный объект к точке.
     * @param LocationInterface|LocationAggregateInterface $location
     * @return LocationInterface
     */
    protected static function convertLocationAggregate( LocationInterface|LocationAggregateInterface $location ):LocationInterface
    {
        if(is_a($location,LocationAggregateInterface::class))
            return $location->getLocation();
        return $location;
    }

    /**
     * Приводит объект LocationInterface|LocationAggregateInterface
     * к массиву вида [longitude, latitude], преобразованных в радианы.
     * @param LocationInterface|LocationAggregateInterface $location
     * @return array{lat:float,long:float}
     */
    protected static function deg2radLocationToArray( LocationInterface|LocationAggregateInterface $location):array
    {
        $location = self::convertLocationAggregate( $location );
        return [
            "lat" => deg2rad( $location->getLatitude() ),
            "long" => deg2rad( $location->getLongitude() )
        ];
    }
}