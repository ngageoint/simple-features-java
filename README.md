# Simple Features Java

#### Simple Features Lib ####

The Simple Features Libraries were developed at the [National Geospatial-Intelligence Agency (NGA)](http://www.nga.mil/) in collaboration with [BIT Systems](http://www.bit-sys.com/). The government has "unlimited rights" and is releasing this software to increase the impact of government investments by providing developers with the opportunity to take things in new directions. The software use, modification, and distribution rights are stipulated within the [MIT license](http://choosealicense.com/licenses/mit/).

### Pull Requests ###
If you'd like to contribute to this project, please make a pull request. We'll review the pull request and discuss the changes. All pull request contributions to this project will be released under the MIT license.

Software source code previously released under an open source license and then modified by NGA staff is considered a "joint work" (see 17 USC § 101); it is partially copyrighted, partially public domain, and as a whole is protected by the copyrights of the non-government authors and must be released according to the terms of the original open source license.

### About ###

[Simple Features](http://ngageoint.github.io/simple-features-java/) is a Java library of geometry objects and utilities based upon the [OGC Simple Feature Access](http://www.opengeospatial.org/standards/sfa) standard.

### Simple Feature Conversion Libraries ###

* [simple-features-wkb-java](https://github.com/ngageoint/simple-features-wkb-java) - Well Known Binary
* [simple-features-geojson-java](https://github.com/ngageoint/simple-features-geojson-java) - GeoJSON
* [simple-features-proj-java](https://github.com/ngageoint/simple-features-proj-java) - Projection

### Usage ###

View the latest [Javadoc](http://ngageoint.github.io/simple-features-java/docs/api/)

### Installation ###

Pull from the [Maven Central Repository](http://search.maven.org/#artifactdetails|mil.nga|sf|2.0.1|jar) (JAR, POM, Source, Javadoc)

```xml

<dependency>
    <groupId>mil.nga</groupId>
    <artifactId>sf</artifactId>
    <version>2.0.1</version>
</dependency>

```

### Build ###

Build this repository using Eclipse and/or Maven:

    mvn clean install
