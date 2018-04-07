CDF   0   
      time             site_id       SGP    facility_id       C1     sds_mode      b1     missing_value         -9999      
sample_int        0.12 seconds   averaging_int         
30 minutes     instruments       0Anemometer:Fill Windmaster Pro, IRGA: LICOR 7500   time_stamp_details        RThe time stamp currently reflects the _start_ of the 30 minute sampling interval.      boom_direction        K0 degrees; wdir is accurate to within an offset of approximately 10 degrees    sign_convention      fc_corr, h, le, mean_g_soil are positive away from surface, r_net is positive towards the surface. That is r_net = h+le+g. Also, r_net = (r_down_short_hemisp+r_down_long_hemisp) - (r_up_short_hemisp+r_up_long_hemisp). fc_corr is defind such that (fc_corr-fc_wpl_h-fc_wpl_le) = fc_uncorr     contact_name      %Marc Fischer, LBNL, mlfischer@lbl.gov      PLEASE_READ_THIS     �For all publications and presentations, please acknowledge: 'U.S. Department of Energy as part of the Atmospheric Radiation Measurement Program.' The automatic inclusion of a data originator as a co-author is not insisted upon in the ARM program, but the source of any data should be clearly recognized either as a co-author or through an appropriate acknowledgment.  The ARM/LBNL Carbon Project contact for this data set is Marc Fischer (mlfischer@lbl.gov). Also please note that we may make adjustments to the data to incorporate adjustments to calibration scales or other issues.  Users should contact Marc Fischer (mlfischer@lbl.gov) to inquire about planned data releases and revisions. Whenever possible, we would appreciate receiving preprints for publications that use the data to insure that the quality and limitations of the data are accurately represented. Your questions and comments are welcome.      MENTOR_QC_FIELD_INFORMATION      �The qc flags use binary, or bitmapped, values (ie values of 0,1,2,4,8,etc.) to note each relevant condition. These can be summed for more detailed qc information. The qc values are:
QC VALUE:  EXPLANATION:
 0-         value not suspect
 1-         value missing
 2-         value below minimum or above maximum or value is +-infinity
 4-         one or more dependencies failed; see the 'dependency' attribute of variable; if a dependency is a 
            spike count, it fails if it is > 100, if a dependency is a variable, it fails when the variable's qc flag
            is not equal to 0
 8-         value has large variance; see 'large variance condition' attribute of variable
16-         value susbect because it has more than 100 spikes; see 'dependency' attribute of variable;
            currently this flag applies to t,q,c
32-         value suspect because of another condition, see 'special condition' attribute of variable.
            Currently only applies to fc_corr, ustar.   	arm_field               irga_serial       C     sonic_serial      B     nrlite_serial         B     history       �created by code 4msonicb1tob1met.c, version r11, with operating system RedHat Linux, kernel 2.4.18-18.7.x, i686 on Apr 09 2015, 08:54:35 GMT      k   	base_time                string        27-Oct-2012,00:00:00 GMT       	long_name         Base time in Epoch     units         #seconds since 1970-1-1 0:00:00 0:00         M�   time_offset                 	long_name         Time offset from base_time     units         'seconds since 27-Oct-2012,00:00:00 GMT          M�   yyyydddhhmmss                   units         yyyydddhhmmss      	long_name         start of integration interval           M�   doy                 units         fractional days    	long_name         fractional day of the year          M�   fc_corr                 units         umol m-2 s-1   	long_name         WPL corrected CO2 flux     	valid_min         ����   	valid_max               
dependency        !nspk_unrot_w,mean_t,mean_q,mean_c           M�   
qc_fc_corr                  units         unitless   	long_name         qc flag for fc_corr         M�   fc_wpl_h                units         umol m-2 s-1   	long_name         #additive WPL H correction to c flux    	valid_min         ����   	valid_max                    M�   qc_fc_wpl_h                 units         unitless   	long_name         qc flag for fc_wpl_h        M�   	fc_wpl_le                   units         umol m-2 s-1   	long_name         $additive WPL LE correction to c flux   	valid_min         ����   	valid_max                    N    qc_fc_wpl_le                units         umol m-2 s-1   	long_name         $additive WPL LE correction to c flux        N   h                   units         W m-2      	long_name         corrected sensible heat flux   	valid_min         ���8   	valid_max               
dependency        nspk_unrot_w,mean_t         N   qc_h                units         unitless   	long_name         qc flag for h           N   le                  units         W m-2      	long_name         WPL corrected latent heat flux     	valid_min         ����   	valid_max               
dependency        nspk_unrot_w,mean_t         N   qc_le                   units         unitless   	long_name         qc flag for le          N   
mean_rot_u                  units         m s-1      	long_name         mean horizontal wind speed          N   qc_mean_rot_u                   units         unitless   	long_name         qc flag for mean_r          N   mean_t                  units         degree C   	long_name         5mean sonic temperature (t), i.e. virtual temperature       	valid_min         ����   	valid_max            2   large_variance_condition      abs(mean_t)/sqrt(var_t) < 2         N    	qc_mean_t                   units         unitless   	long_name         qc flag for mean_t          N$   mean_q                  units         mmol m-3   	long_name         mean H2O density   	valid_min         <#�
   	valid_max         D���   large_variance_condition      abs(mean_q)/sqrt(var_q) < 2         N(   	qc_mean_q                   units         unitless   	long_name         qc flag for mean_q          N,   mean_c                  units         mmol m-3   	long_name         mean CO2 density   	valid_min         A (�   	valid_max         A��        N0   	qc_mean_c                   units         unitless   	long_name         qc flag for mean_c          N4   mean_p                  units         kPa    	long_name         mean IRGA pressure     	valid_min            ^   	valid_max            e        N8   	qc_mean_p                   units         unitless   	long_name         qc flag for mean_p          N<   	var_rot_u                   units         (m s-1)2   	long_name         variance of u      	valid_min                	valid_max            d        N@   qc_var_rot_u                units         unitless   	long_name         qc flag for var_r           ND   	var_rot_v                   units         (m s-1)2   	long_name         variance of v      	valid_min                	valid_max            d        NH   qc_var_rot_v                units         unitless   	long_name         qc flag for var_r           NL   	var_rot_w                   units         (m s-1)2   	long_name         variance of w      	valid_min                	valid_max                    NP   qc_var_rot_w                units         unitless   	long_name         qc flag for var_r           NT   var_t                   units         C2     	long_name         variance of t      	valid_min                	valid_max                    NX   qc_var_t                units         unitless   	long_name         qc flag for var_t           N\   var_q                   units         (mmol m-3)2    	long_name         variance of q      	valid_min                	valid_max           p        N`   qc_var_q                units         unitless   	long_name         qc flag for var_q           Nd   var_c                   units         (umol m-3)2    	long_name         variance of c      	valid_min                	valid_max         >���        Nh   qc_var_c                units         unitless   	long_name         qc flag for var_c           Nl   wdir                units         degrees    	long_name         horizontal wind direction      	valid_min                	valid_max           h        Np   qc_wdir                 units         unitless   	long_name         qc flag for wdir        Nt   theta                   units         degrees    	long_name         rotation to zero w     	valid_min                	valid_max           h        Nx   qc_theta                units         unitless   	long_name         qc flag for theta           N|   phi                 units         degrees    	long_name         rotation to zero v     	valid_min                	valid_max           h        N�   qc_phi                  units         unitless   	long_name         qc flag for phi         N�   ustar                   units         m s-1      	long_name         friction velocity      	valid_min                	valid_max            
        N�   qc_ustar                units         unitless   	long_name         qc flag for ustar           N�   Lmoni                   units         meters     	long_name         Monin-Obukhov length scale     	valid_min         ����   	valid_max           #(        N�   qc_Lmoni                units         unitless   	long_name         qc flag for Lmoni           N�   bar_pres                units         kPa    	long_name         barometric pressure    	valid_min            ^   	valid_max            e        N�   qc_bar_pres                 units         unitless   	long_name         qc flag for bar_pres        N�   t_air_upper                 units         degree C   	long_name         !mean air temperature, upper level      	valid_min         ����   	valid_max            2        N�   qc_t_air_upper                  units         unitless   	long_name         qc flag for t_air_upper         N�   t_air_lower                 units         degree C   	long_name         !mean air temperature, lower level      	valid_min         ����   	valid_max            2        N�   qc_t_air_lower                  units         unitless   	long_name         qc flag for t_air_lower         N�   rh_upper                units         %      	long_name         #mean relative humidity, upper level    	valid_min                	valid_max            d        N�   qc_rh_upper                 units         unitless   	long_name         qc flag for rh_upper        N�   rh_lower                units         %      	long_name         #mean relative humidity, lower level    	valid_min                	valid_max            d        N�   qc_rh_lower                 units         unitless   	long_name         qc flag for rh_lower            N�   z_upper              units         meters     	long_name         height of upper RH/T sensor         M�   z_lower              units         meters     	long_name         height of lower RH/T sensor         M�   mean_g_soil                 units         W m-2      	long_name         -average of up to four soil heat flux sensors       	valid_min         ����   	valid_max            �   special_condition         3less than two out of four sensors read valid values         N�   qc_mean_g_soil                  units         unitless   	long_name         qc flag on mean_g_soil          N�   stderr_mean_g_soil                  units         none   	long_name         standard error of mean_g_soil           N�   z_g_soil             units         meters     	long_name         depth of soil heat flux probes          M�   r_net                   units         W m-2      	long_name         net solar radiation (NR-lite)      	valid_min         ����   	valid_max                    N�   qc_r_net                units         unitless   	long_name         qc flag for r_net           N�   r_tot                   units         W m-2      	long_name         5total downwelling shortwave radiation (LI-200 sensor)      	valid_min         ����   	valid_max           �        N�   qc_r_tot                units         unitless   	long_name         qc flag for r_tot           N�   r_down_short_hemisp                 units         W m-2      	long_name         -downwelling shortwave hemispheric irradiance       	valid_min         ����   	valid_max           �        N�   qc_r_down_short_hemisp                  units         unitless   	long_name         qc flag for r_down_short_hemisp         N�   r_up_short_hemisp                   units         W m-2      	long_name         *upwelling shortwave hemispheric irradiance     	valid_min         ����   	valid_max           �        N�   qc_r_up_short_hemisp                units         unitless   	long_name         qc flag for r_up_short_hemisp           N�   r_down_long_hemisp                  units         W m-2      	long_name         <downwelling longwave hemispheric irradiance; t_rad corrected   	valid_min         ����   	valid_max               
dependency        t_rad           N�   qc_r_down_long_hemisp                   units         unitless   	long_name         qc flag for r_down_long_hemisp          N�   r_up_long_hemisp                units         W m-2      	long_name         :upwelling longwave hemispheric irradiance; t_rad corrected     	valid_min         ����   	valid_max               
dependency        t_rad           N�   qc_r_up_long_hemisp                 units         unitless   	long_name         qc flag for r_up_long_hemisp        N�   t_rad                   units         degree K   	long_name         5radiometer body temperature (for longwave correction)      	valid_min            �   	valid_max           J   
dependency        t_rad           N�   qc_t_rad                units         unitless   	long_name         qc flag for t_rad           O    ppfd                units         umol m-2 s-1   	long_name         (photosynthetic photon flux density (PAR)   	valid_min         ����   	valid_max           �        O   qc_ppfd                 units         unitless   	long_name         qc flag for ppfd        O   mean_t_soil_upper                   units         degree C   	long_name         -average of two upper soil temperature sensors      	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O   qc_mean_t_soil_upper                units         unitless   	long_name         qc flag for mean_t_soil_upper           O   stderr_mean_t_soil_upper                units         none   	long_name         Sstandard error of mean_t_soil_upper; note maximum of two values used in calculation         O   mean_t_soil_middle                  units         degree C   	long_name         .average of two middle soil temperature sensors     	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O   qc_mean_t_soil_middle                   units         unitless   	long_name         qc flag for mean_t_soil_middle          O   stderr_mean_t_soil_middle                   units         none   	long_name         Tstandard error of mean_t_soil_middle; note maximum of two values used in calculation        O    mean_t_soil_lower                   units         degree C   	long_name         -average of two lower soil temperature sensors      	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O$   qc_mean_t_soil_lower                units         unitless   	long_name         qc flag for mean_t_soil_lower           O(   stderr_mean_t_soil_lower                units         none   	long_name         Sstandard error of mean_t_soil_lower; note maximum of two values used in calculation         O,   mean_m_soil_upper                   units         cm3 cm-3   	long_name         +average of four upper soil moisture sensors    	valid_min                	valid_max         ?      special_condition         /less than two of four sensors read valid values         O0   qc_mean_m_soil_upper                units         unitless   	long_name         qc flag for mean_m_soil_upper           O4   stderr_mean_m_soil_upper                units         none   	long_name         $standard error of mean_m_soil_upper         O8   mean_m_soil_lower                   units         cm3 cm-3   	long_name         +average of four lower soil moisture sensors    	valid_min                	valid_max         ?      special_condition         /less than two of four sensors read valid values         O<   qc_mean_m_soil_lower                units         unitless   	long_name         qc flag for mean_m_soil_lower           O@   stderr_mean_m_soil_lower                units         none   	long_name         $standard error of mean_m_soil_lower         OD   precip                  units         
mm 30min-1     	long_name         mean precipitation     	valid_min            ^   	valid_max            e        OH   	qc_precip                   units         unitless   	long_name         qc flag for precip          OL   lat              units         degrees    	long_name         latitude of instrument location         M�   lon              units         degrees    	long_name          longitude of instrument location        M�   alt              units         meters     	long_name         #height of tower base from sea level         M�   ppfd_up                 units         umol m-2 s-1   	long_name         2Upwelling photosynthetic photon flux density (PAR)     	valid_min         ����   	valid_max           �        OP   
qc_ppfd_up                  units         unitless   	long_name         qc flag for ppfd_up         OT   zm               units         	meters         	long_name         #height of instrument from tower bas         M�   zrad             units         meters     	long_name         height of radiation instruments         M�   z_t_soil_upper               units         meters     	long_name         depth of upper t soil probes        M�   z_t_soil_middle              units         meters     	long_name         depth of middle t soil probes           M�   z_t_soil_lower               units         meters     	long_name         depth of lower t soil probes        M�   z_m_soil_upper               units         meters     	long_name         depth of upper m soil probes        M�   z_m_soil_lower               units         meters     	long_name         depth of lower m soil probes        M�P�$ @|(�@	���L��Bk����jC�� @���@   �L�;����  ���;L��        B}Hg��  C�� ��ӫ  ��W��  �9�MW  ����  ��*��  �@L��  �Aza�  �C�V  �AV�  �B�+�  �>_�  �=�N�  �=98  �=ص�  �>o�  �:�0  �Av��  �>�$  �C$�%  �>%��  �A��L  �B��  ��ff  �@�S�  �B~�  �B2
=  ��S��    @n�Rff  ��<  ���b  �?���  �C��}  �C��:  �C���  ���-  �@��  �=��@���  �>"�@�~�  �;�� �<  !��< �<  !��<       ����  �@�      B}Hgȏ� C���<��  ���sQ  ��&�d  ��ܞ9  ��q)  �@,!k  �An�	  �C�K:  �A�y�  �B�/B  �='�  �<��   �;�do  �=g�	  �?�_  �9�.�  �C��6  �=�^  �C1�i  �=r`t  �@-�P  �B�)  ���I�  �@��  �B�  �BFff  ���    @t#��3  ��<  �����  �?��  �C�n'  �C���  �C��f  ��49X  �@t�  �=�dP?���  �>G�@��  �=�n��<  !��< �<  !��<       ��-��  �@�      B}Hg�E  C��U���  ��"I�  ��S�  ��^��  ����  �@7s  �Al�  �C�  �A�Ŵ  �B�1@  �>B  �=+  �<��(  �=��  �@S.  �9�'�  �C���  �?D۲  �C+EL  �=L�  �>��E  �B�$Z  ��	�  �@u�-  �B�    �BL�
  ����P    @���zQ�  ��<  ����  �@{o  �C��  �C���  �C�&f  ��9X  �?��;  �>n�?�X  �>�I�?��  �=����<  !��< �<  !��<       ��t�  �@�     B}Hg� � C�� :���  ��S)  ����  ���x  ��;�  �@x  �Af�&  �C��  �A�t5  �B�1�  �>3�  �=�;�  �<�=�  �>&:$  �@��  �;W�  �C�Ɏ  �?��  �C�  �=�
}  �@ȇ;  �B�(�  ��&�y  �@U�  �B�33  �BV�  ���P    @���m�  ��<  ����F  �@h��  �C}[  �C��  �C�ٚ  ��Ƨ�  �@ (�  �>!G�?�hs  �>��@ bN  �; �<  !��< �<  !��<       �����  �@�      B}Hg̶  C����$�&  ���Wc  �:�{  ��	��  �>��  �?�|�  �A[��  �C�Z  �A�TE  �B�1�  �>��  �<� }  �:��  �=�T�  �A�q  �:zA�  �C�Y  �?�&�  �C�  �<��  �?���  �B�/  ��cdZ  �@ �`  �B���  �Bc�  ����T    @q$�bff  ��<  ����h  �?ޗ�  �C}'  �C���  �C�Y�  �����  �>��  �>��h>š�  �>�33?G�  �=u�<  !��< �<  !��<       ���hs  �@��     B}Hg�q� C��U��j�  ��v��  �9�Wv  ���/�  ��~�  �?�Z%  �A^0&  �C���  �A�>.  �B�6�  �<�Y�  �;�&�  �:�O�  �= ��  �?�Q�  �9�&*  �C��t  �?��t  �C`�  �<�  �?v�  �B�1'  ��a%  �@  �B���  �Bi  �����    @\�[�P�\  ��<  ��!x�  �@���  �Czr�  �C��S  �C��   ���  �����  �?|��h�9  �?��:��9X  �=N��<  !��< �<  !��<       ��6E�  �@�     B}Hg�'  C�� ��K�  ��hR  �;{��  ��LJ!  �>�M�  �?�z  �ART}  �C�_�  �A���  �B�42  �=S�L  �<ml&  �:��  �=�Sa  �@���  �<Ҿ  �C��  �?���  �C�  �=_�  �@�aw  �B�/  ����u  �?��  �B�ff  �Br�\  ��'\,    @�s��Y��  ��<  ���+  �@A��  �Cy�^  �C��  �C��3  ����P  ��b  �>��C�GK�  �?녿�  �=�n��<  !��< �<  !��<       ���j  �@Ȝ     B}Hg�� C������[  ����  ��P_�  ��
�  ���s  �?��  �AL��  �C���  �A�\g  �B�1n  �<�Z�  �<k  �;��  �=�JG  �@сd  �;�@  �C��  �?��$  �C�<  �=�   �?�  �B�)y  ����  �>��7  �B���  �B}�H  ���^6    @rS[�]Q�  ��<  ��?��  �@O�;  �Cw�[  �C���  �C�@   ��H�9  �����  �>�+���  �>�ȴ��V  �=� ��<  !��< �<  !��<       �=�{  �@�      B}Hgј  C��U�n_  ��WD�  �9��+  ����  ��*�  �?�Kr  �AL   �C��Z  �A�Q  �B�1+  �<�y�  �<��  �:f��  �=nƲ  �?�=�  �;,�  �C���  �?Hy  �C��  �<�1R  �?��  �B�%�  ���j  �>�  �B�    �B��R  ��n�    @�#��T33  ��<  ��5�h  �@~�+  �Cvy�  �C�I�  �C�&f  ��o  ����`  �>�Ĝ�ۥ�  �>�����n  �=m�p�<  !��< �<  !��<       �=�F  �@Ϥ     B}Hg�S� C�� ���k  ��efL  �:.��  ���8k  �=��  �?�Q  �AA�O  �C�u�  �A��  �B�1  �=��  �=!�6  �:+og  �=�  �A  �;��  �C��V  �?1�  �C��  �<��  �@wK  �B�%�  ����  ��&$�  �B�    �B���  ��;b    @c���M��  ��<  ����  �@8 �  �Ct  �C���  �C��f  �;�o  ��,�  �>V�(�2ff  �>�^4�*ff  �=J���<  !��< �<  !��<       �>��  �@є     B}Hg�	  C���?V�  ��	�  ����  ��6�  ��C�H  �@)�m  �A=�0  �C���  �A���  �B�,!  �=�=  �=)"�  �<44N  �>�  �A[�  �;^8  �C��  �?�^  �C	�  �=�=  �@�W�  �B��  ���    ��V  �B�33  �B���  ��=�    @�j�[��  ��<  ��J-  �@g|�  �Cr?�  �C���  �C��   �<���  ����R  �>���l  �>~�ؿ���  �<����<  !��< �<  !��<       �>��!  �@�V     B}Hg�Ā C��U>O�T  ��Z�  ��#��  ���v�  ���	�  �@"��  �A74v  �C�*�  �A��  �B�$  �=n��  �=��2  �<Q��  �>	�q  �A1l�  �;g�  �C�B�  �?���  �C�  �=�_s  �@"*x  �B�+  ���"�  ����P  �B�    �B���  ��d��    @]���]Q�  ��<  ����  �@�-  �Cr�  �C�9  �C��   �<�j  ���^  �=�M����  �>-X�  �<�@�<  !��< �<  !��<       �>��  �@�     B}Hg�z  C�� ��yl  ����  �;P��  ����  �>���  �?��  �A-g7  �C��  �A��^  �B�&�  �=�'C  �=�g�  �;\!�  �=�p�  �A	'�  �;:N  �C���  �?�KF  �B��	  �<�.�  �?)^�  �B�  ���p�  ��	�#  �B���  �B���  ��`��    @P&��\    ��<  �����  �?�M�  �Cq��  �C���  �C�Y�  �=H�9  ��U��  �>(���W�  �>E���C�w  �>����<  !��< �<  !��<       �?ff  �@��     B}Hg�5� C������`  ���;u  ���Y�  ���$�  ��s��  �?�  �A,��  �C��>  �A�\�  �B�&D  �=���  �=nwf  �<�  �=�E�  �@M)  �;�  �C�=f  �?:�  �Cz�  �=M3)  �@�7  �B�  ���\)  ���
  �B�    �B���  ��A��    @\���T�  ��<  ��V5?  �@E�  �Cp��  �C���  �C��  �>J  ��h�j  �>G+�l�t  �>�M��_�<  �>#��<  !��< �<  !��<       �?3t�  �@؜     B}Hg��  C��U����  ���Ω  �90��  ��	�  �<�a  �?��R  �A*��  �C�]A  �A�G�  �B�)  �=7��  �<��4  �:���  �<�f�  �@[��  �:ME  �C�#�  �?�{  �B��V  �<�\�  �B�<  �B��  ���ff  ��'�w  �B���  �B�ff  ���bM    @K�X�
  ��<  ����F  �@O
=  �Cl��  �C�1P  �C��f  �>hr�  ���hs  �>b����bN  �>MO���l�  �=_;��<  !��< �<  !��<       �?L1  �@�^     B}Hg٦� C�� ��   ��Y�7  �;#]I  ����8  �>��{  �@Ț  �A.�  �C��P  �A��x  �B�'b  �=7�  �>2�f  �;��  �=�`~  �@�Tb  �:�ct  �C�]�  �?2U:  �C
�  �=�%  �E�"�  �B��  ���  ��r�  �B�    �B�33  ��Ȳ    @i���P33  ��<  ��	�^  �@�S�  �Cn��  �C��S  �C��f  �>�V  ��n~�  �>)x��v$�  �>�=t�n  �=\ �<  !��< �<  !��<       �?i7L  �@�      B}Hg�\  C�����$�  ��"�:  �;r�  ��Yb  �>��w  �?�ES  �A3C�  �C�pT  �A�Gc  �B�  �<�@�  �=�(  �;��  �<�z�  �@|"�  �:9�  �@>�   �?��  �C1�  �=X{p  �@} �  �B��?  �����  �����  �Bș� �B�ff  ��.��    @)'�_ff  ��<  ��F�+  �@��F  �Co\a  �C��U  �C��3  �?�  ����/  �>�C�����  �?�x���  �=!���<  !��< �<  !��<       �?kƨ  �@��     B}Hg�� C��U�j9  �=:;  �;�c   �?+�|  �?.�  �?j
� �A+cc  �C��  �A��	  �B�   �=��  �>o��  �;6��  �=�x�  �A#�|  �;  �C�q�  �>���  �C�  �<L�  ��{��  �B��  ���
=  ��n�  �B�ff �B�    ��333    @���R��  ��<  ��,I�  �@m?}  �Cm*�  �C�g�  �C��3  �?%�  ���{  �>V� ���  �>�hp���  �<�\ �<  !��< �<  !��<       �?�`B  �@ߤ     B}Hg��  C�� =���  ���Z  ��I�c  ��2L|  ���e~  �?���  �A&��  �C��y  �A�^�  �B���  �=��o  �=��  �<eT  �=��f  �@�-b  �:���  �C�Cz  �>��W  �C��  �=���  �@f�S  �B��Z  �����  ��H��  �B�   �B���  �����    @%��H33  ��<  ��M�  �@�J  �Cm��  �C�¾  �C�s3  �?��9  ���?}  �>�0��2  �>�����  �<�`��<  !��< �<  !��<       �?�  �@�     B}Hgވ� C������  �� Z  �<6�J  ��S:�  �?�#*  �?�M2  �A��  �C��a  �A��x  �B���  �=�z>  �>t�{  �</�  �>��|  �B�P  �<&��  �C�g�  �?�  �C�*  �=o�X  �A�~  �B��  ��  ��E�-  �Bԙ� �B���  ��z�    @/�Y�Vz�  ��<  ��6$�  �@ff  �Cl�T  �C��9  �C��f  �?���  ���p�  �>V�0��X  �>�j����H  �;Ĝ �<  !��< �<  !��<       �?��  �@�     B}Hg�>  C��U>@�   ��-��  �;�  ��f7�  �=���  �?���  �A#��  �C��F  �A�Lu  �B���  �=$��  �>7�  �;��  �=�(P  �@��  �;X�  �C��|  �>���  �C+�  �<��	  �>�7�  �B��V  ���
  ��Y7L  �B�   �B�    ��[��    @���W33  ��<  ��	G�  �@|j  �ClH  �C�	O  �C��  �?�v�  ���/  �>������`  �?&ȴ���  �=����<  !��< �<  !��<       �?�+  �@�u     B}Hg��� C�� ����  �<�c  �;�KP  �>�"�  �?3�  �?��!  �A��  �C���  �A�B�  �B���  �=��  �="Z�  �:A�  �=���  �AN��  �;�y  �C��n  �>�χ  �B�i�  �=+�\  ��KO  �B��  ���  ���X  �Bݙ� �B�    ��+S�    @1���K33  ��<  ��bM�  �@��  �Ck  �C�0�  �C��3  �?��-  ���J  �>��@��Z  �>�����9X  �=�- �<  !��< �<  !��<       �?�j  �@�V     B}Hg�  C���=�,  ��O��  �:��  ���̄  �=�  �?��  �A�F  �C��  �A��K  �B���  �<�;  �;��  �:�/  �=�q  �@H�%  �;�=  �C�)P  �?Cg  �B�@  �<ة�  �@Ϭ�  �B��B  ��
=  ����  �B��� �B�    ����8    @[$��J=q  ��<  ���\  �@t�j  �Cj��  �C���  �C���  �?��y  ��ʴ9  �>�`@���y  �>�����l�  �=�n��<  !��< �<  !��<       �?���  �@�7     B}Hg�j� C��U���y  �=��i  �;�1  �?���  �?(�  �?⥲  �A ��  �C��  �A���  �B�:  �<�N�  �;Ȟ�  �:E�  �=��5  �@�C8  �:�<�  �C��  �?`�,  �B��:  �<�Z�  ��7T�  �B���  ��    ���S�  �B♚ �B���  ��<�B    @J�N�I�
  ��<  ��x�  �@:n�  �Ci��  �C���  �C��3  �?��-  ��У�  �>S���Гu  �>k� ��&�  �={�@�<  !��< �<  !��<       �?�l�  �@�     B}Hg�   C�� ��i�  �����  �:���  ���n�  �>-`�  �?�Q�  �@��)  �C���  �A��+  �B� 4  �=1l�  �<��\  �:��  �=J�~  �@� �  �:v*�  �C�+  �?@��  �B�4�  �<ɫ  �?|��  �B���  ���H  �����  �B�ff �B���  ��p�    @J���H    ��<  ��!&�  �@q&�  �CiH�  �C�@�  �C��f  �?�+  ��ӡ�  �>"M���n�  �>BP����  �=.��<  !��< �<  !��<       �?�  �@��     B}Hg�ۀ C�«���~  ��x��  ����V  ����  ��F�@  �@ g  �@���  �C���  �A��_  �B�#g  �=�ԡ  �<�O  �;�p&  �=g��  �@�\�  �:2/�  �C�ͤ  �?Oғ  �C�G  �=��  �?ʻ  �B��  ���H  ����u  �B�   �B���  ���    @��>�R  ��<  �@f�R  �@�v�  �Cj�B  �C�^�  �C���  �A;
=  ���{  �>�`�Ɨ�  �>49P�Å  �=\@�<  !��< �<  !��<       �@Cƨ  �@��     B}Hg�  C��U���t  ��
��  ����  ��32�  ���+  �?�;  �@���  �C��[  �A��x  �B�)�  �>=�  �<�d�  �;��  �=��z  �A��  �:�P\  �C��  �?"�b  �B�O�  �=Jiz  �@v��  �B��  ���\  �����  �B�ff �B���  ��_    A����z�  ��<  �BO�  �A�(�  �Ci��  �C�Da  �C��3  �B�    ��E�^  �>�j��m�h  �>��<�J-  �=�S��<  !��< �<  !��<       �A��  �@�     B}Hg�L� C�� =��  �>�1[  �<���  �A��  �@P��  �?�S   �A��  �C���  �A��  �B�/�  �>c�_  �=��>  �<���  �>�S  �B.��  �;`~  �C���  �>�w�  �C��  �=�v�  �� c�  �B��  ��\)  ��49X  �B�   �B�    �B���    A^fdA�  ��<  �C��  �B    �Cl?>  �C��Z  �C�    �CV    �>V  �?G����  �>�r���+  �>�u�<  !��< �<  !��<       �B+�H  �@�     B}Hg�  C�ʫ?O�  �?�5_  �=hֳ  �B�=  �@�j�  �@�k  �A�  �C��  �A�˩  �B�/�  �>"�}  �?(  �={n  �>B��  �BB_  �;Z`t  �C�2\  �>���  �Cd�  �>�F  ����  �B�'�  ����H  ��j  �Bϙ� �B�33  �B��    Aϙ�B|{  ��<  �Cy�   �B]�H  �Cp�r  �C��l  �C�s3  �C��f  �@9��  �?
~�@0bN  �?
~�@j  �=\@�<  !��< �<  !��<       �B�Ǯ  �@�}     B}Hg꽀 C��U@H,  �@693  �=�v�  �Br��  �A_��  �@��  �A+Z  �C��3  �A��e  �B�+C  �?�V  �>�;�  �=���  �>�q;  �A���  �<F�0  �A���  �>���  �C#�  �>b"#  ��|�  �B�/  ����w  �@C�  �B���  �B��R  �µ��B   �< C��  ��<  �C��  �B�ff  �Cw�  �C��  �C��  �D ��  �@�=p  �?Q�@���  �>��P@���  �?M��<  !��< �<  !��<       �B�    �@�^     B}Hg�s  C�� ?�l  �@��k  �>
��  �B�u�  �A��j  �@T �  �AAY  �C�h  �A���  �B��  �?�  �?�W�  �=��0  �>�wn  �@�;d  �;KV  �B;�0  �      �C�  �>r�n  ��&I  �B��  ��$j  �@x1'  �B���  �Br\)  ��< B  �< CL�f  ��<  �C�&f  �B���  �C~�  �C���  �C�ff  �DN    �@���  �?~5D@�\)  �>��@�=p  �>���<  !��< �<  !��<       �B���  �@�?     B}Hg�.� C�ҫ@'�  �@�  �=���  �C ��  �A���  �@�!�  �AO��  �C�V�  �A��6  �B��  �?:��  �?o��  �>'$�  �>��J  �@�(�  �;\�/  �B�G�  �>�65  �Bٸ  �>���  ��65^  �B��j  ���`B  �@�ƨ  �B���  �Ba��  ��< B  �< C��  ��<  �D    �B���  �C��T  �C���  �C�33  �Dv�   �Az�  �?uAz�  �?%�A �H  �?%�;�s  ��< �<  !��<       �C33  �@�      B}Hg��  C��U@A�L  �A&Y�  �>�  �Cd��  �A�Xv  �@VH�  �AU�  �C�c�  �A�~  �B��h  �?S�  �?�>  �>#�  �?'t�  �A]  �;8V>  �B��D  �?��}  �B�E�  �>��  ���x�  �B��?  ��p��  �@�(�  �B�33  �BP�  ��< B  �< C��3  ��<  �D��  �B�    �C�E  �C��  �C��f  �D��   �A
�H  �=#� A��  �?Q�Az�  �?nx<+�N  ��< �<  !��<       �C	��  �@�     B}Hg C�� ? �  �A1f�  �>'/�  �Ct��  �A���  �@77� �AU��  �C��L  �A�3  �B��	  �?�Y  �?��  �>-��  �?I�2  �@�є  �;ߊ  �B���  �?��  �B�Yb  �>sTU  ����  �B��  ���A�  �@�I�  �B���  �BEp�  �CF33B   �< C���  ��<  �D*`   �B�    �C��{  �C��  �C���  �D�@   �A�H  �<�� A)  �?�PA5�  �=�@<�G  ��< �<  !��<       �CL�  �@��     B}Hg�U  C�ګ>Xv&  �A��  �=�r�  �C@mO  �A��9  �@.�u �AYa�  �C� C  �A�E�  �Bř  �?��f  �?l�~  �>�  �?8�  �A��A  �;��  �B�ZH  �?~T=  �B���  �>�!O  ��q��  �BŷL  �>49X  �@�Q�  �B��  �B9�  �C(��B   �< Cǳ3  ��<  �D2@   �B�33  �C� �  �C��>  �C�    �D�@   �A:(�  �?XQ�A5  �>�(�A5  �>��P<��|  ��< �<  !��<       �C�3  �@��     B}Hg�� C��U@���  �A"��  �>o�)  �C`��  �B�K  �@Hg �A_��  �C���  �A�O  �BŇ�  �?��P  �@��  �>4n�  �?H�#  �B  �;O<�  �B�B�  �@%'  �B��  �>���  ��muk  �BŘ�  �?;dZ  �@�
=  �B�B�  �B/�  �C33B   �< CӀ   ��<  �D8�   �B�    �C��  �C�8�  �C�s3  �D��   �AI
>  �?��AH��  �?#�A?�
  �=�� <���  ��< �<  !��<       �C�   �@�     B}Hg��  C�� @W߳  �A6?  �=�s�  �C}\  �A�i;  �@� �Af�  �C�.�  �A~+/  �B�k�  �?Ϸ�  �?�f  �>��  �?w�]  �A�1  �;?�  �Bċ�  �?l��  �B�t  �>=�U  ��l�  �B�}q  �?���  �A
�R  �BlQ�  �B  �B�ffB   �< C�L�  ��<  �D;�   �B虚  �C���  �C�R�  �C�ٚ  �D�`   �AS
>  �?�G�AO�  �?�fhA@��  �?.x<�^J  ��< �<  !��<       �C��  �@�B�    B}Hg� C��=���  �A�'�  �=�*�  �C���  �A���  �@	�� �Ao��  �C�2�  �A}M�  �B�9p  �?�,  �?��  �>M6A  �?��8  �A �  �;Z/N  �B�t  �@9i  �B���  �>�*f  ���O�  �B�Y  �@
�\  �A�  �B^z�  �B��  �B�p    ?� Cԙ�  ��<  �D9�   �B�33  �C�  �C�b  �C�@   �D�    �AY��  �?ǮAX(�  �?��PAQ�  �?@  <���  ��< �<  !��<       �C33  �@�     B}Hg�7  C��U�s��  �A%հ  �>��  �Cg�P  �AҞD  �?��l �Ar�>  �C���  �A|�  �B�.�  �?�&�  �@Bk  �>kD  �?@�e  �A��s  �;�&  �B���  �@,�4  �BZF'  �>��  ����  �B�=�  �@��  �AG�  �BU\)  �Bz�  �@��`    A^oC���  ��<  �D3�   �B���  �C��  �C�V  �C�s3  �D��   �A`(�  �?��HAm�  �?��A[�  �?y���<  !��< �<  !��<       �B�ff  �@�#�    B}Hg�� C�� ���  �A�D  �>��m  �CP�  �B%Nn  �?�m� �Aw�!  �C��p  �A|rM  �B�z  �?�n�  �@(m�  �>�K  �?B  �A�Uh  �;��  �B��h  �@C  �BV�2  �>��  �����  �B�/�  �@A&�  �A&ff  �BDG�  �B	=q  ��ծ    A�SaC��   ��<  �D$�   �B�ff  �C���  �CϜ�  �C���  �D��   �A[��  �?�34Ai34  �?���A\Q�  �?ffh<��  ��< �<  !��<       �B�    �@�     B}Hg��  C���2�=  �@ۋ*  �=��`  �C�
  �Ai��  �?�� �Az�"  �C��  �A{�>  �B���  �?�-�  �>�"E  �>'c  �?V�  �@���  �;��  �C��  �A�a  �B�g  �?N�  �¬�-  �B��  �@lz�  �A2�R  �B/z�  �B (�  ��
>    A��	C�L�  ��<  �D��  �B���  �C��A  �C��"  �C�ٚ  �D��   �Ah=p  �?�Q�A|�\  �@)G�A~�\  �=�\ <�"S  ��< �<  !��<       �B�    �@��    B}Hg�c� C��U�M,�  �@���  �=�b�  �CG�  �Ak�2  �?zA �A�#�  �C��  �A{��  �B��h  �?^��  �?@�  �=�{�  �>�b  �@��j  �:�&�  �B�*T  �@�ˏ  �B�լ  �>��  ��NO  �B��q  �@�M�  �A9G�  �B�
  �A�G�  �����    A��Cd�3  ��<  �C�33  �B�    �C���  �C���  �C��3  �DX@   �Ab��  �?���A��  �@��Az=p  �>�
0<���  ��< �<  !��<       �B�    �@�u     B}Hg�  C�� >��  �@��  ��Ny  �B���  �@l-�  �?J� �A�=�  �C��  �A{^�  �B�ע  �?�х  �?Gw�  �=��  �>��I  �A�B  �:l�  �B1|�  �@-��  �C��  �>���  ���p  �B��  �@��  �A6�R  �B!�  �A�=q  ���yX    AlQmC9L�  ��<  �CҌ�  �B���  �C�sJ  �C�r  �C��f  �D;�   �Aj��  �?H��A�z�  �?���Ay��  �>P<��M  ��< �<  !��<       �B�    �@��    B}Hg�Ԁ C��5��  �?��  �=V��  �B(��  �@���  �?t�� �A�f  �C���  �A{6{  �B���  �>�e�  �?���  �=���  �>T��  �@V�H  �9�'w  �B�H  �@��  �B��  �=b��  ����  �B���  �@w�  �A3
=  �B$33  �A��  ���j    Au�B�aH  ��<  �Cq33  �B!  �C��  �C��  �C�@   �C��  �AS�  �>��Aep�  �?n�A[\)  �>W
@<��Z  ��< �<  !��<       �B0Q�  �@�V     B}Hg��  C��U�&�6  ��N�  �=r�  ���l�  �@{D  �?"4 �A��j  �C��  �A{�  �B��  �>Р6  �>���  �=JCS  �>��  �?��  �9��  �B�O�  �@̤m  �B{`�  �>tRW  �C���  �B��`  �@��D  �A7�  �B�  �A�{  ��#C�    A%�B,z�  ��<  �C[L�  �B,�  �C��  �C�\T  �C��   �C��f  �A\Q�  �<�� At    �?���Ah�R  �>=p�<���  ��< �<  !��<       �B*�  �@�ƀ    B}Hg�E� C�� �>�  �>ܜ+  �;���  �A,  �?uA  �?��� �A�Z+  �C��P  �A{*{  �B�×  �>��  �>���  �=�  �=��   �? ��  �9>�%  �B��x  �@Fۣ  �B��  �>>�  ��l�!  �B�ݲ  �@y��  �A5G�  �B"\)  �A�R  ���    A��r�'|�  ��<  �Cff  �A���  �C���  �C��  �C��f  �Cw��  �AKp�  �>ٙ�Adz�  �?�fhAX��  �>8R <�`l  ��< �<  !��<       �A�33  �@�7     B}Hg��  C����+  ���  ���  ��=U>  ���/G  �?��4  �A�-;  �C��B  �A{��  �B���  �=��i  �=ō  �;i�  �=���  �@�D�  �9�  �BʉH  �>̖b  �B�v�  �<��  �>��T  �B��  �@O�w  �A#
=  �B'(�  �A�{  ����    @�A��F{  ��<  �B2�  �A1G�  �C��  �C��  �C���  �B�    �A�  �>�f`A(��  �>�f`A �H  �=��<��  ��< �<  !��<       �A33  �@���    B}Hh �� C��U<�Ҝ  ����4  �:*�  ���S  �;�4  �?w� �A���  �C��>  �A|q�  �B���  �=lI�  �>���  �:���  �=?��  �?�  �:{w|  �C^  �?���  �A�^  �<���  �>��  �B�ٚ  �@	��  �A�H  �B?��  �B�  ���Ĝ    ?�NI�m
=  ��<  ��t�  �?��  �C���  �C�/9  �C���  �@�    �@��:  �?%`D@��  �?��@��`  �<�C�<,��  ��< �<  !��<       �?
��  �