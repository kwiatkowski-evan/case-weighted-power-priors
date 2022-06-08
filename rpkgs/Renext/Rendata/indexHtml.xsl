<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0">
  <!-- =============================================================
       Author: Y. Deville
       Style sheet to display a description of the 'index.xml' file
       in a HTML table
       
       Simple modifications could produce LaTeX code and thus pdf.

    ================================================================ -->
  <xsl:output method="html" encoding="iso-8859-1"/>
  <xsl:template match="/">
    <html>
      <body style="margin-left:100;margin-right:1000;font-size:8pt">
        <h1>
	  Indexed content of the directory
	</h1>
	<p>
	  Each dataset contains <strong>one</strong><code> OTdata</code> and optionally
	  one or several historical data <code>OTSdata</code> or
	  <code>MAXdata</code>. While <strong>OTdata</strong> is the main 
	  piece of information, both historical data types are designed
	  to contain extra information.    
	</p>
	<ul>
	  <li>
            A <code>MAXdata</code> contains the largest observations
	    for events in a specified period. This corresponds to the 
	    "<i>r</i>-largest" data type. The number of events
	    must be <i>r</i>&#160;&gt;=1.
          </li>
	  <li>
            A <code>OTSdata</code> contains the observations
	    for events with level above a specified threshold which must be
	    larger than the <code>OTdata</code> threshold. It can contain
	    zero events (<i>r</i>&#160;=0), in which case  the provided information 
	    is that the OTS threshold
	    was never exceeded during the specified period.
          </li>
	  <li>
	    <code>OTdata</code> may contain <strong>missing periods</strong>, but 
	    neither <code>OTSdata</code> nor <code>MAXdata</code>
	    can.
	  </li>
	  <li>
	    <code>OTdata</code> may contain events or missing periods
	    read from a <strong>csv file</strong>, but 
	    neither <code>OTSdata</code> nor <code>MAXdata</code>
	    can.
	  </li>
	</ul>
	<p>
	  Datatime format can be chosen for data read from files, 
	  but should within the XML index conform to 
	  "YYYY-mm-dd HH:MM:SS".
        </p>
        <table border="1" width="880" cellspacing="0"
               style="font-size:8pt border-collapse: collapse"
               bordercolor="#C0C0C0">
          <xsl:attribute name="cellpadding">4</xsl:attribute>
          <xsl:attribute name="bgcolor">#F7F7F7</xsl:attribute>
	  <tr>
	    <th>name<br/>shortlab</th>
	    <th>var</th>
	    <th>(main)<br/>OTdata</th>
	    <th>(historical)<br/>OTSdata</th>
	    <th>(historical)<br/>MAXdata</th>
	    <th>describe</th>
	  </tr>
          <xsl:for-each select="datasets/dataset">  
            <tr>
              <td valign="top">
		<strong><xsl:value-of select="@name"/></strong><br/>
		<xsl:value-of select="@shortLab"/>
              </td>
              <td valign="top">
		<strong><xsl:value-of select="@varName"/></strong><br/>
		(<xsl:value-of select="@varUnit"/>)
              </td>
              <td valign="top" width="280"> 
		<xsl:for-each select="OTdata">
		  <font color="SteelBlue"><b>info</b></font><br/>
		  <font color="grey">start </font> 
		  <xsl:value-of select="@start"/><br/>
		  <font color="grey">end </font> 
		  <xsl:value-of select="@end"/><br/>
		  <font color="grey">effDuration </font>
		  <xsl:value-of select="@effDuration"/><br/>
		  <font color="grey">threshold </font>
		  <xsl:value-of select="@threshold"/><br/><br/>
		  <font color="SteelBlue"><b>data</b></font><br/>
	          <xsl:for-each select="data">
		    <xsl:for-each select="file">
		      <font color="grey">file </font>
		      <a href="{@path}">
			<xsl:value-of select="@path"/>
		      </a><br/>
		    </xsl:for-each>
		    <xsl:for-each select="events">
		      <font color="grey">events</font>
		      (<xsl:value-of select="count(event)"/>)<br/>
		      <xsl:for-each select="event[position() &lt;= 2]">
			  <xsl:value-of select="@date"/><br/>
		      </xsl:for-each>
		      <xsl:choose>
			<xsl:when test="count(event) &gt; 2">
			  ...
			</xsl:when>
		      </xsl:choose>
		    </xsl:for-each>
                  </xsl:for-each> <br/>
		  <font color="SteelBlue"><b>missing</b> </font><br/>
		  <xsl:for-each select="missing">
		    <xsl:for-each select="file">
		      <font color="grey">file </font>
		      <a href="{@path}">
			<xsl:value-of select="@path"/>
		      </a><br/>
		    </xsl:for-each>
		    <xsl:for-each select="periods">
		      <font color="grey">periods </font>
		      (<xsl:value-of select="count(period)"/>)<br/>
		      <xsl:for-each select="period[position() &lt;= 2]">
			<xsl:value-of select="@start"/> 
			<font color="grey"> to </font>
			<xsl:value-of select="@end"/><br/>
		      </xsl:for-each>
		      <xsl:choose>
			<xsl:when test="count(period) &gt; 2">
			  ...
			</xsl:when>
		      </xsl:choose>
		    </xsl:for-each>
		  </xsl:for-each>
		</xsl:for-each>
	      </td>
	      <td valign="top" width="180"> 
		<font color="SteelBlue"><b>OTSdata</b> </font>
		(<xsl:value-of select="count(OTSdata)"/>)<br/>
		<xsl:for-each select="OTSdata[position() &lt;= 2]">
		  <font color="grey">shortLab </font> 
		  <xsl:value-of select="@shortLab"/><br/>
		  <font color="grey">start </font> 
		  <xsl:value-of select="@start"/><br/>
		  <font color="grey">end </font> 
		  <xsl:value-of select="@end"/><br/>
		  <!-- <font color="grey">Dur. </font>
		  <xsl:value-of select="@duration"/><br/> -->
		  <font color="grey">threshold </font>
		  <xsl:value-of select="@threshold"/><br/>
		  <font color="grey">evts r </font> 
		  <xsl:value-of select="count(data/events/event)"/><br/><br/>
		</xsl:for-each>
		<xsl:choose>
		  <xsl:when test="count(OTSdata) &gt; 2">
		    ...
		  </xsl:when>
		</xsl:choose>
              </td>
              <td valign="top" width="180"> 
		<font color="SteelBlue"><b>MAXdata</b> </font>
		(<xsl:value-of select="count(MAXdata)"/>)<br/>
		<xsl:for-each select="MAXdata[position() &lt;= 2]">
		  <font color="grey">shortLab </font> 
		  <xsl:value-of select="@shortLab"/><br/>
		  <font color="grey">start </font> 
		  <xsl:value-of select="@start"/><br/>
		  <font color="grey">end </font> 
		  <xsl:value-of select="@end"/><br/>
		  <!-- <font color="grey">Dur. </font>
		  <xsl:value-of select="@duration"/><br/> -->
		  <font color="grey">evts r </font> 
		  <xsl:value-of select="count(data/events/event)"/><br/><br/>
		</xsl:for-each>
		<xsl:choose>
		  <xsl:when test="count(MAXdata) &gt; 2">
		    ...
		  </xsl:when>
		</xsl:choose>
              </td>
              <td valign="top" width="380">
		<xsl:copy-of select="describe"/>
              </td>
            </tr>
          </xsl:for-each>
	</table>
      </body>
    </html>
  </xsl:template>
</xsl:stylesheet>





