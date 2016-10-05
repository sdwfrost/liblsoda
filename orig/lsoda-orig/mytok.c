char *mytok(s1, s2)
register unsigned char *s1;
register unsigned char *s2;
{
	static unsigned char *cp1;

	if (s1 != NULL)
		cp1 = s1;
	s1 = cp1 + strspn(cp1, s2);
	if (*s1 == '\0')
		return(NULL);
	cp1 = s1 + strcspn(s1, s2);
	*cp1++ = '\0';
	return(s1);
}
