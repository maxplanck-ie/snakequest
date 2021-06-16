FROM centos:7.2.1511

MAINTAINER Thomas Manke, manke@ie-freiburg.mpg.de

RUN yum install -y epel-release && \
    yum install -y wget ypbind yp-tools ypserv autofs nfs-utils rsyslog && \
    yum install -y openssl-devel curl libcurl-devel mesa-libGLU libpng-devel cairo-devel php && \
    yum install -y R-3.6.0-1.el7.x86_64 && \
    yum clean all && \
    R -e 'install.packages(c("shiny","shinydashboard","shinyBS","rhandsontable","reshape2","sendmailR"), repos="https://cran.rstudio.com/",dependencies=TRUE, clean=TRUE)'

RUN R -e 'install.packages(c("shinyalert"), repos="https://cran.rstudio.com/",dependencies=TRUE, clean=TRUE)'

# Fix the time zone
RUN unlink /etc/localtime && \
    ln -s /usr/share/zoneinfo/Europe/Berlin /etc/localtime

RUN mkdir /data /root/snakequest /etc/automount

COPY mounts.py /usr/local/bin/mounts.py
COPY startup.sh /usr/local/bin/startup.sh
COPY app.R /root/snakequest/app.R
COPY Rprofile.site /usr/lib/R/etc/Rprofile.site

VOLUME ["/export/"]

EXPOSE 2525

CMD ["/usr/local/bin/startup.sh"]



