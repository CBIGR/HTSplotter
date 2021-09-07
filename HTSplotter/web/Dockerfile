FROM php:8.0.3-apache

LABEL "maintainer"="jasper.anckaert@ugent.be"
LABEL "description"="HTSPlotter"
LABEL "version"="0.1"

# Update OS and install common dev tools
RUN apt-get update
RUN apt-get install -y wget nano vim git zip unzip zlib1g-dev libzip-dev libpng-dev libjpeg-dev libhdf5-serial-dev

# Install python3 and modules
RUN apt-get install -y python3 python3-pip
RUN pip3 install numpy
RUN pip3 install matplotlib
RUN pip3 install scipy
RUN pip3 install h5py
RUN pip3 install minio

# Install PHP extensions needed
RUN docker-php-ext-install -j$(nproc) gd zip pcntl exif

# Enable common Apache modules
RUN a2enmod headers expires rewrite

# Install Composer
COPY --from=composer:2 /usr/bin/composer /usr/local/bin/composer

# Set user to www-data
USER www-data

# Create tmp_dir
RUN mkdir /var/www/html/user_files
RUN chmod 766 /var/www/html/user_files

#Add php WebUI code
COPY ./web/ /var/www/html/

#Add python application scripts
RUN mkdir /var/www/html/scripts
COPY *.py /var/www/html/scripts/

# Set working directory to workspace
WORKDIR /var/www/html/
