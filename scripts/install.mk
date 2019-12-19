exec_prefix:=$(prefix)

define POST_INSTALL_MSG
	@echo "--------- Library installed in $(prefix)."
	@echo "* Add the following line to '$(HOME)/.bashrc' to use pkg-config"
	@echo "* PKG_CONFIG_PATH=${prefix}/lib/pkgconfig:PKG_CONFIG_PATH"
	@echo "* The command make post-install will make that for you"
	@echo "* --------------------------------------------------------------"
	@echo "* Then you can use pkg-config to compile your program"
	@echo "* See \"Use instructions\" in README.md for more information"
endef

pcfile:=$(PRJ)_$(compiler).pc
INSTALL_LIB= $(exec_prefix)/lib
INSTALL_INCLUDE= $(prefix)/include/$(PRJ)
INSTALL_PKG=$(INSTALL_LIB)/pkgconfig

install: $(numlib) $(SCRPTD)/$(pcfile)
	@echo "-------------------------------------------------------------"
	@echo "--------- Installing in $(prefix)"
	@echo "-------------------------------------------------------------"
	$(MKDIR) $(INSTALL_LIB) $(INSTALL_INCLUDE) $(INSTALL_PKG)
	$(CP) $(numlib) $(INSTALL_LIB)/
	$(CP) $(MODD)/$(PRJ).mod $(INSTALL_INCLUDE)/
	$(CP) $(SCRPTD)/$(pcfile) $(INSTALL_PKG)/
	cd $(INSTALL_PKG); $(LN) $(pcfile) $(PRJ).pc; cd -
	@echo "-------------------------------------------------------------"
	$(call POST_INSTALL_MSG)

uninstall:
	$(RM) -r $(INSTALL_INCLUDE)
	$(RM) $(INSTALL_LIB)/lib$(PRJ).a
	$(RM) $(INSTALL_PKG)/$(PRJ).pc $(INSTALL_PKG)/$(pcfile)
	@echo "* Remember that you may have to update your $$HOME/.bashrc file"

# pkg-config file for specific compiler	
$(SCRPTD)/$(pcfile): $(SCRPTD)/$(PRJ).pc.in
	sed -e 's|@PRJ@|${PRJ}|' -e 's|@version@|${VERSION}|' -e 's|@prefix@|${prefix}|' -e 's|@FC@|${FC}|' $< > $@


post-install:
	@echo 'export PKG_CONFIG_PATH="$(prefix)/lib/pkgconfig:$${PKG_CONFIG_PATH}"' >> $(HOME)/.bashrc
	@source ~/.bashrc
