exec_prefix:=$(prefix)

 POST_INSTALL_MSG
--------- Library installed in $(prefix). 

Add the following line to $HOME/.bashrc to use pkg-config to compile

PKG_CONFIG_PATH=${prefix}/lib/pkgconfig:PKG_CONFIG_PATH"
endef

export POST_INSTALL_MSG

pcfile:=$(PRJ)_$(compiler).pc
INSTALL_LIB= $(exec_prefix)/lib
INSTALL_INCLUDE= $(prefix)/include/$(PRJ)
INSTALL_PKG=$(prefix)/pkgconfig
install: $(numlib) $(SCRPTD)/$(pcfile)
	@echo "-------------------------------------------------------------"
	@echo "--------- Installing in $(prefix)"
	@echo "-------------------------------------------------------------"
	$(MKDIR) $(INSTALL_LIB) $(INSTALL_INCLUDE) $(INSTALL_PKG)
	$(CP) $(numlib) $(INSTALL_LIB)/
	$(CP) $(MODD)/$(PRJ).mod $(INSTALL_INCLUDE)/
	$(CP) $(SCRPTD)/$(pcfile) $(INSTALL_PKG)/
	$(LN) $(INSTALL_PKG)/$(pcfile) $(INSTALL_PKG)/$(PRJ).pc
	@echo "-------------------------------------------------------------"
	@echo "$$POST_INSTALL_MSG"


# pkg-config file for specific compiler	
$(SCRPTD)/$(pcfile): $(SCRPTD)/$(PRJ).pc.in
	sed -e 's|@PRJ@|${PRJ}|' -e 's|@version@|${VERSION}|' -e 's|@prefix@|${prefix}|' -e 's|@FC@|${FC}|' $< > $@




post-install:
	echo "PKG_CONFIG_PATH=${prefix}/lib/pkgconfig:PKG_CONFIG_PATH"
